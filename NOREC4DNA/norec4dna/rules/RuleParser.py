#!/usr/bin/python
# -*- coding: latin-1 -*-
import typing
from re import compile, search
from collections import Counter

from norec4dna.helper.fallback_code import strContainsSub_python, microsatellite_python, longestSequenceOfChar_python

try:
    import cdnarules
except ModuleNotFoundError:
    print("C Module failed to load, falling back to slow mode")


def microsatellite(text: typing.AnyStr, length_to_look_for: int) -> typing.Tuple[int, str]:
    try:
        return cdnarules.microsatellite(text, length_to_look_for)
    except NameError:
        return microsatellite_python(text, length_to_look_for)


def longestSequenceOfChar(text: typing.AnyStr, char_x="*") -> typing.Tuple[str, int]:
    try:
        return cdnarules.longestSequenceOfChar(text, char_x)
    except NameError:
        return longestSequenceOfChar_python(text, char_x)


def strContainsSub(text: typing.AnyStr, sequence: typing.AnyStr) -> bool:
    try:
        return cdnarules.strContainsSub(text, sequence)
    except NameError:
        return strContainsSub_python(text, sequence)


debug = False

"""
data  = string (combination of A, C, G, T) e.g. AAGTCAGAAGTGTTGAAAAAAAAAGCGTGTTTGCC
rules = List of Rules, a Rule is a Tuple of ("ruleKind(parameter)" (as String), "dropProbability" (in Range (0.0,1.0] ))
"""

""" ALL RULES MUST HAVE THE (temporary) FORM rule(x,y,z) """


# @jit
def switch(name: str) -> typing.Callable[[typing.AnyStr, typing.Any, typing.Any], int]:
    switcher = {
        "longestSequenceOfChar": (
            lambda x, y, z: 1 if int(z) <= longestSequenceOfChar(x, y)[1] else 0
        ),
        "strContainsSub": (lambda x, y, z: 1 if strContainsSub(x, y) else 0),
        "strContainsSubRegex": (
            lambda x, y, z: 1 if strContainsSubRegex(x, y) else 0
        ),
        "strContainsIllegalChars": (
            lambda x, y, z: 1 if strContainsIllegalChars(x, y) else 0
        ),
        "charCountBiggerEqualThanX": (
            lambda x, y, z: 1 if int(z) <= charCountBiggerEqualThanX(x, y) else 0
        ),
        "microsatelliteLongerThanX": (
            lambda x, y, z: 1 if int(z) <= microsatellite(x, int(y))[0] else 0
        ),
        "len": (lambda x, y, z: 1 if int(y) <= length(x) else 0),
        "gcContent": (lambda x, y, z: 1 if float(y) >= gc_content(x) else 0),
        "gcContentLQ": (lambda x, y, z: 1 if float(y) <= gc_content(x) else 0),
        "*": (lambda x, y, z: 1),
    }
    return switcher.get(name)


def gc_content(text: typing.AnyStr) -> float:
    counter = Counter(text)
    count = counter["G"] + counter["C"]
    return (count / len(text)) * 100


def iupac_replace(sequence: typing.AnyStr):
    iupac_regex = {'M': '[AC]', 'R': '[AG]',
                   'W': '[AT]', 'S': '[CG]',
                   'Y': '[CT]', 'K': '[GT]',
                   'V': '[ACG]', 'H': '[ACT]',
                   'D': '[AGT]', 'B': '[CGT]',
                   'X': '[ACGT]', 'N': '[ACGT]'}
    for i, j in iupac_regex.items():
        sequence = sequence.replace(i, j)
    if debug:
        print(sequence)
    return compile(sequence)


def strContainsSubRegex(text: typing.AnyStr, sequence: typing.AnyStr) -> bool:
    iupac_seq = iupac_replace(sequence)
    res = search(iupac_seq, text)
    if debug:
        print(res)
    return bool(res)


def strContainsIllegalChars(text: typing.AnyStr, allowed_chars: typing.AnyStr) -> int:
    for cha in text:
        if cha not in allowed_chars:
            return 1
    return 0


def charCountBiggerEqualThanX(text: typing.AnyStr, cha: typing.AnyStr):
    res = text.count(cha)
    if debug:
        print(res)
    return res


def length(text: typing.AnyStr) -> int:
    res = len(text)
    if debug:
        print(res)
    return res


def executeRule(rule_kind: str, data: typing.AnyStr) -> int:
    if "(" not in rule_kind:
        rule_kind += "(*,0)"
    name, params = rule_kind.split("(")
    params = params.split(")")[0].replace(" ", "")
    if params == "":
        params = "*,0"
    if "," not in params:
        params += ",0"
    p1, p2 = params.split(",")
    if p2 == "":
        p2 = 0
    return switch(name)(data, p1, p2)


def shouldDrop(data: typing.AnyStr, rules: typing.List[typing.Tuple[str, float]]) -> float:
    drop_chance = 0.0
    for rule in rules:
        rule_kind, drop_prob = rule
        drop_chance += executeRule(rule_kind, data) * drop_prob
    return min(1.0, drop_chance)


def shouldDropMax(data: typing.AnyStr, rules: typing.List[typing.Tuple[str, float]]) -> float:
    drop_chance = 0.0
    for rule in rules:
        rule_kind, drop_prob = rule
        tmp = executeRule(rule_kind, data) * drop_prob
        if drop_chance < tmp:
            drop_chance = tmp
    return min(1.0, drop_chance)


def shouldDropMin(data: typing.AnyStr, rules: typing.List[typing.Tuple[str, float]]) -> float:
    drop_chance = 1.0
    for rule in rules:
        rule_kind, drop_prob = rule
        tmp = executeRule(rule_kind, data) * drop_prob
        if drop_chance > tmp:
            drop_chance = tmp
    return min(1.0, drop_chance)


if __name__ == "__main__":
    print(microsatellite("ACGACGAAGAAGAAGAAGAGTAGTAGAAGA", 2))

    in_text = "AAAGCCGAGAGAATTTTTTCACAAAAAAAAAAAAAGTTATAAATCCAATCA" * 10000
    rules1 = [
        ("*", 0.1),
        ("len(15)", 0.2),
        ("strContainsIllegalChars(ACGT)", 1.0),
        ("charCountBiggerEqualThanX(A, 15)", 0.5),
        ("microsatelliteLongerThanX(2,10)", 0.01),
    ]
    print(
        shouldDrop(in_text, rules1)
        + shouldDropMax(
            in_text,
            [
                ("longestSequenceOfChar(*,2)", 0.001),
                ("longestSequenceOfChar(*,3)", 0.005),
                ("longestSequenceOfChar(*,4)", 0.01),
                ("longestSequenceOfChar(*,5)", 0.05),
                ("longestSequenceOfChar(*,6)", 0.1),
                ("longestSequenceOfChar(*,7)", 0.9),
                ("longestSequenceOfChar(*,8)", 1.0),
            ],
        )
    )
    print(
        shouldDrop(
            in_text,
            [
                ("gcContent(50)", 0.001),
                ("gcContent(40)", 0.01),
                ("gcContent(30)", 0.02),
                ("gcContent(20)", 0.03),
                ("gcContent(10)", 0.04),
                ("gcContent(0)", 0.05),
            ],
        )
    )
    print(
        shouldDrop(
            in_text,
            [
                ("strContainsSub(TATAAA)", 0.01),
                ("strContainsSub(TTGACA)", 0.05),
            ],
        )
    )
    print(
        shouldDrop(
            in_text,
            [
                ("strContainsSubRegex(CANYYY)", 0.01),
                ("strContainsSubRegex(ANCCAATCA)", 0.01),
            ],
        )
    )
