"""SYSTEM.md Hard Rules loading."""

from agents.common import system_rules, with_system_rules


def test_system_rules_content():
    rules = system_rules()
    assert "Hard Rules" in rules
    assert "pseudobulk" in rules.lower()
    assert "300 dpi" in rules.lower()
    assert "provenance" in rules.lower()


def test_with_system_rules_prepends():
    wrapped = with_system_rules("agent prompt body")
    assert wrapped.startswith("# scAgent Hard Rules")
    assert "agent prompt body" in wrapped
