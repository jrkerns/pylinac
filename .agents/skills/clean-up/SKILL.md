---
name: clean-up
description: Clean up code with lightweight clarity improvements
disable-model-invocation: true
---

# clean-up

When cleaning up, the goal is to make simple modifications that improve clarity
or structure of the code without modifying lots of logic.

Don't be aggressive
Opportunities for refactor include:
* Adding or clarifying types
* Adding or clarifying docstrings
* Adding comments to code blocks
* Using newer syntax supported by the language
* Clarifying documentation

Do not make helper functions that are 2 lines or less. In that case, don't refactor.
