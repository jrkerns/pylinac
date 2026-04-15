# AGENTS.md: Pylinac User Docs Guide

## Purpose

Source-level guide for authored user documentation in `docs/source/`.

This file applies only to `docs/source/**`. Ignore generated output such as `docs/build/**`, `docs/build-pdf/**`, and any other build artifacts that may appear after local docs builds.

Refer to `../../AGENTS.md` for repo-wide conventions and development workflow.

## Source Intent Map

| Path | Intent |
| --- | --- |
| `docs/source/*.rst` | Top-level user documentation pages and module guides |
| `docs/source/topics/` | Cross-cutting conceptual and how-to documentation |
| `docs/source/images/` | Referenced images, HTML artifacts, and media used by pages |
| `docs/source/files/` | Referenced sample data files used in examples or docs assets |
| `docs/source/code_snippets/` | Python snippets included with `literalinclude` |
| `docs/source/conf.py` | Source of truth for enabled Sphinx extensions and docs behavior |

## Authoring Rules

- Use ReStructuredText for authored documentation pages. Do not add Markdown docs under `docs/source/`.
- Write for end users first. Assume the reader is a medical physicist who needs practical guidance more than implementation detail.
- Match the existing docs style: clear section headings, explanatory prose, concrete examples, and focused `note` or `warning` callouts where needed.
- Use explicit directive languages such as `.. code-block:: python` and `.. code-block:: bash`.
- Use Sphinx cross-references where helpful, such as `:ref:`, `:class:`, `:func:`, and `:mod:`, instead of loose prose references.
- Add new pages to the appropriate `toctree` when needed. `index.rst` is the top-level navigation entry point.
- Prefer patterns already used in this tree over inventing new presentation styles.
- Do not add `.. contents::` to normal pages unless there is an explicit reason. The current Furo styling does not handle duplicated in-page tables of contents well.

## RST House Style

- When a sentence ends with a colon and is followed by a list, add a blank line before the first list item.
- Keep lists readable and visually separated from the lead-in sentence, even if a tighter form might still render.
- Keep directive indentation exact. Most rendering issues in this tree come from indentation mistakes.

Example:

```rst
Supported options:

- First option
- Second option
```

## Available Syntax

Available docs syntax should come from the configured docs stack, not generic Sphinx assumptions. Check `pyproject.toml` under `[project.optional-dependencies].docs` and `docs/source/conf.py` before using unfamiliar directives.

The currently available stack includes:

- `sphinx` for core Sphinx directives and cross-references
- `furo` for theme styling
- `sphinx-copybutton` for copy-button support on code blocks
- `sphinx-design` for design directives, including `.. tab-set::` and `.. tab-item::`
- `autodoc-pydantic` for Pydantic API documentation support
- `matplotlib.sphinxext.plot_directive` for `.. plot::`

Prefer syntax already present in `docs/source/**`. If a directive is not installed and enabled here, do not assume it is available.

## Validation

After modifying files in `docs/source/**`:

- Run `uvx prek run --files <changed source paths>`.
- Run `uv run nox -s build_docs` for substantive page edits, new pages, new directives, `literalinclude` changes, or `toctree` changes.

Current pre-commit checks relevant to docs include:

- `rst-backticks`
- `rst-directive-colons`
- `rst-inline-touching-normal`
- `blacken-docs`

## Maintenance Protocol

Update this file whenever:

- The `docs/source/` structure changes
- Docs dependencies change in `pyproject.toml`
- Enabled Sphinx extensions change in `docs/source/conf.py`
- New recurring house-style rules should be enforced for authored docs

Keep this guide focused on authored source files only. Do not expand it into a guide for generated docs output.
