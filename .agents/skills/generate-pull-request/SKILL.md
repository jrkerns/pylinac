---
name: generate-pull-request
description: Generate a draft Bitbucket pull request with proper formatting
disable-model-invocation: true
---

# generate-pull-request

Generating a Pull Request is against the pylinac BitBucket repository.
PRs should be created in "draft" mode, not in "ready" mode.
Unless otherwise stated, the PR should be against the "master" branch.

## Git

The branch name should ideally include the JIRA ticket, but is not required. Prompt me
if it does not contain a ticket reference before creating a PR; it's possible there isn't one.

## Code Check

If there is a ticket, it will almost always require a changelog entry in `changelog.rst`.
If the changelog entry is missing, prompt me to add it before creating the PR. In cases
of infrastructure changes such as CI or tooling, no changelog entry is required.

## PR Requirements

- 2 reviewers
- Title: ticket + short summary (e.g., "RAM-4567 Proof of concept of using Python .whl files")
- Description: ticket reference + summary. Screenshots for user-facing features.
- Author comments on lines with non-obvious logic are encouraged.

## References

- Branch naming examples and PR description templates: `references/pr-examples.md`
