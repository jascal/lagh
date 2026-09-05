"""A restricted mathematical grammar for caller-DECLARED forms (jascal/lagh#1).

`sp.sympify(text)` is `eval` in a costume: the string is handed to Python's
parser and evaluated, so any Python expression inside it RUNS -- before a single
mathematical check -- with the permissions of the process that called it.
Measured: the MCP `verify` tool passed its caller-controlled `form` straight in,
and `"(__import__('builtins').print('X'), Symbol('x_0'))[1]"` printed, then
certified `x_0`. Swapping `sympify` for another eval-based parser (`parse_expr`
with transformations, `eval` with a restricted namespace) keeps the hole: the
text still reaches the interpreter.

So the text never reaches the interpreter. It is parsed to a syntax tree
(`ast.parse`, which executes nothing), and the tree is walked EXPLICITLY:
every node kind that is not on the list below is rejected by name, and the
SymPy expression is built node by node from constructors we chose. There is no
namespace to escape from, because nothing is looked up by name except the
declared variables, a handful of named constants, and a fixed function table.

Allowed:
    variables      x_0 .. x_{d-1} (the symbols the caller passes)
    constants      E, pi, GoldenRatio, EulerGamma, Catalan
    numbers        integer and float literals (no complex, no strings)
    operators      + - * / **   (and ^ as power, rewritten to ** before parsing)
    unary          - +
    functions      sqrt cbrt exp log ln sin cos tan asin acos atan atan2
                   sinh cosh tanh asinh acosh atanh Abs abs sign
                   Rational Integer Float (numeric literal arguments only)
    grouping       parentheses

Rejected, always: attribute access, subscripts, calls to anything not in the
table, keyword arguments, tuples/lists/dicts/sets, comparisons, boolean
operators, lambdas, comprehensions, f-strings, walrus, imports, names not in
the tables, and any literal that is not a plain int or float. Resource bounds
(node count, nesting depth, numeric exponent size) keep a pathological form
from being a denial of service instead of a payload.
"""

from __future__ import annotations

import ast

import sympy as sp

MAX_TEXT = 4000            # characters
MAX_NODES = 500            # syntax-tree nodes
MAX_DEPTH = 40             # nesting
MAX_NUMERIC_EXPONENT = 1000  # |exponent| when both base and exponent are numbers


class FormError(ValueError):
    """The declared form is not an expression of the restricted grammar."""


_CONSTANTS = {
    "E": sp.E, "pi": sp.pi, "GoldenRatio": sp.GoldenRatio,
    "EulerGamma": sp.EulerGamma, "Catalan": sp.Catalan,
}

# name -> (constructor, min_arity, max_arity, numeric_args_only)
_FUNCTIONS = {
    "sqrt": (sp.sqrt, 1, 1, False), "cbrt": (sp.cbrt, 1, 1, False),
    "exp": (sp.exp, 1, 1, False), "log": (sp.log, 1, 2, False),
    "ln": (sp.log, 1, 1, False),
    "sin": (sp.sin, 1, 1, False), "cos": (sp.cos, 1, 1, False),
    "tan": (sp.tan, 1, 1, False),
    "asin": (sp.asin, 1, 1, False), "acos": (sp.acos, 1, 1, False),
    "atan": (sp.atan, 1, 1, False), "atan2": (sp.atan2, 2, 2, False),
    "sinh": (sp.sinh, 1, 1, False), "cosh": (sp.cosh, 1, 1, False),
    "tanh": (sp.tanh, 1, 1, False),
    "asinh": (sp.asinh, 1, 1, False), "acosh": (sp.acosh, 1, 1, False),
    "atanh": (sp.atanh, 1, 1, False),
    "Abs": (sp.Abs, 1, 1, False), "abs": (sp.Abs, 1, 1, False),
    "sign": (sp.sign, 1, 1, False),
    "Rational": (sp.Rational, 1, 2, True), "Integer": (sp.Integer, 1, 1, True),
    "Float": (sp.Float, 1, 1, True),
}

_BINOPS = {
    ast.Add: lambda a, b: a + b,
    ast.Sub: lambda a, b: a - b,
    ast.Mult: lambda a, b: a * b,
    ast.Div: lambda a, b: a / b,
    ast.Pow: None,          # handled explicitly (exponent bound)
}


def _reject(node, what: str | None = None):
    raise FormError(f"{what or type(node).__name__} is not allowed in a declared form")


class _Builder:
    def __init__(self, syms):
        self.names = {str(s): s for s in syms}
        self.count = 0

    def build(self, node, depth: int = 0):
        self.count += 1
        if self.count > MAX_NODES:
            raise FormError(f"form has more than {MAX_NODES} nodes")
        if depth > MAX_DEPTH:
            raise FormError(f"form is nested deeper than {MAX_DEPTH}")
        if isinstance(node, ast.Expression):
            return self.build(node.body, depth + 1)
        if isinstance(node, ast.Constant):
            return self._constant(node)
        if isinstance(node, ast.Name):
            return self._name(node)
        if isinstance(node, ast.UnaryOp):
            v = self.build(node.operand, depth + 1)
            if isinstance(node.op, ast.USub):
                return -v
            if isinstance(node.op, ast.UAdd):
                return v
            _reject(node.op, f"unary {type(node.op).__name__}")
        if isinstance(node, ast.BinOp):
            return self._binop(node, depth)
        if isinstance(node, ast.Call):
            return self._call(node, depth)
        _reject(node)

    def _constant(self, node):
        v = node.value
        # bool is an int subclass and is NOT a number here
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            _reject(node, f"literal of type {type(v).__name__}")
        try:
            return sp.Integer(v) if isinstance(v, int) else sp.Float(v)
        except Exception as e:                                 # noqa: BLE001
            raise FormError(f"bad numeric literal {v!r}: {e}") from e

    def _name(self, node):
        if node.id in self.names:
            return self.names[node.id]
        if node.id in _CONSTANTS:
            return _CONSTANTS[node.id]
        raise FormError(
            f"unknown symbol {node.id!r}: variables are "
            f"{', '.join(self.names) or 'x_0..x_{d-1}'} and constants "
            f"{', '.join(_CONSTANTS)}")

    def _binop(self, node, depth):
        op = type(node.op)
        if op not in _BINOPS:
            _reject(node.op, f"operator {op.__name__}")
        a = self.build(node.left, depth + 1)
        b = self.build(node.right, depth + 1)
        try:
            if op is ast.Pow:
                if isinstance(b, sp.Number) and isinstance(a, sp.Number) \
                        and abs(float(b)) > MAX_NUMERIC_EXPONENT:
                    raise FormError(
                        f"numeric exponent {b} exceeds {MAX_NUMERIC_EXPONENT}")
                return a ** b
            return _BINOPS[op](a, b)
        except FormError:
            raise
        except Exception as e:                                 # noqa: BLE001
            raise FormError(f"cannot build {op.__name__}: {e}") from e

    def _call(self, node, depth):
        if not isinstance(node.func, ast.Name):
            _reject(node.func, "calling anything but a named function")
        name = node.func.id
        if name not in _FUNCTIONS:
            raise FormError(
                f"unknown function {name!r}: allowed are {', '.join(_FUNCTIONS)}")
        if node.keywords:
            _reject(node, "keyword arguments")
        fn, lo, hi, numeric_only = _FUNCTIONS[name]
        if not (lo <= len(node.args) <= hi):
            raise FormError(f"{name} takes {lo}-{hi} argument(s), got {len(node.args)}")
        args = [self.build(a, depth + 1) for a in node.args]
        if numeric_only and not all(isinstance(a, sp.Number) for a in args):
            raise FormError(f"{name} takes numeric literals only")
        try:
            return fn(*args)
        except Exception as e:                                 # noqa: BLE001
            raise FormError(f"cannot build {name}: {e}") from e


def parse_form(text: str, syms) -> sp.Expr:
    """Parse `text` as an expression of the restricted grammar over `syms`.

    Returns a SymPy expression built node by node, or raises `FormError`.
    Nothing in `text` is ever evaluated as Python."""
    if not isinstance(text, str):
        raise FormError(f"form must be a string, got {type(text).__name__}")
    if len(text) > MAX_TEXT:
        raise FormError(f"form longer than {MAX_TEXT} characters")
    # `^` means power, as sympify's convert_xor did -- rewritten at the TEXT
    # level, as that transformation was, so it binds like `**` (as a BitXor
    # node it would bind looser than `*`: x^2*y would parse as x**(2*y))
    text = text.strip().replace("^", "**")
    try:
        tree = ast.parse(text, mode="eval")
    except (SyntaxError, ValueError) as e:
        raise FormError(f"not a mathematical expression: {getattr(e, 'msg', e)}") from e
    expr = _Builder(syms).build(tree)
    if not isinstance(expr, sp.Basic):
        raise FormError("form did not build a mathematical expression")
    return sp.sympify(expr)   # already a SymPy object: this is the identity, no text
