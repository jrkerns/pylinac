# PR Review Checklist

## Bitbucket PR Metadata

* The PR title should reference a ticket if applicable. Example: "RAM-1234 Improve docs".
* The PR description should list the problem being solved or the feature being added. Referencing or copying from the ticket is also acceptable if the ticket contains that information.
* The Bitbucket pipelines should be passing.
* PRs should have at least 2 reviewers assigned under normal circumstances.
* The task list should be complete.
* Take into consideration the PR description. These can be valuable to rationalize design choices.
* Take into consideration any Bitbucket PR comments. These can be valuable to rationalize design choices.
* The branch should be up to date with master.

## PR Scope

* PRs should be focused on a single change, bug fix, or feature. Flag unrelated fixes for separate PRs.
* Target less than 300 lines of added/modified code (net change), excluding tests and documentation.

## Code Quality

* Review the algorithms for reasonable correctness and performance.
* Review the code for readability and maintainability.

## Algorithms & Performance

* If an algorithm or expression is inefficient or could be improved, provide actionable alternatives. Example:
  - Bad: "This is too slow."
  - Good: "I think this can be done better by doing X instead of Y. What do you think?"

## Tests

* There should always be tests for customer-facing features and bug fixes. If there are no tests it should be explained in the PR description why there aren't.
* Tests should cover the "happy path" as well as edge cases and be as comprehensive as possible.

## Documentation

* Customer-facing features and bugfixes should include a note in `changelog.rst`.

## Security

* Check for any security vulnerabilities that might be introduced by this code.
* No private or sensitive information should be committed to the repository.
