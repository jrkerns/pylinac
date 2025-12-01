# pull-request-review

## General

* If this is a BitBucket PR, feel free to pull down the PR branch locally for better context.
* If this is a Bitbucket PR, be sure to examine any BitBucket PR comments.
* Do not add the review as a comment in the PR. Let me review it locally before posting.

## Bitbucket PR metadata

* The PR title should reference a ticket if applicable in the PR title and/or description. Example: "RAM-1234 Improve docs".
* The PR description should list the problem being solved or the feature being added. Referencing or copying from the ticket is also acceptable if the ticket contains that information.
* The bitbucket pipelines should be passing.
* PRs should have at least 2 reviewers assigned under normal circumstances.
* The task list should be complete.
* Take into consideration the PR description. These can be valuable to rationalize design choices.
* Take into consideration any Bitbucket PR comments. These can be valuable to rationalize design choices.
* The branch should be up to date with master.

## Code Review

### Size, scope

* Review the algorithms for reasonable correctness and performance.
* Review the code for readability and maintainability.
* Beautiful is better than ugly.
  Explicit is better than implicit.
  Simple is better than complex.
  Complex is better than complicated.
  Flat is better than nested.
  Sparse is better than dense.
  Readability counts.
  Special cases aren't special enough to break the rules.
  Although practicality beats purity.
  Errors should never pass silently.
  Unless explicitly silenced.
  In the face of ambiguity, refuse the temptation to guess.
  There should be one-- and preferably only one --obvious way to do it.
  Although that way may not be obvious at first unless you're Dutch.
  Now is better than never.
  Although never is often better than *right* now.
  If the implementation is hard to explain, it's a bad idea.
  If the implementation is easy to explain, it may be a good idea.
  Namespaces are one honking great idea -- let's do more of those!*
* PRs should be as small as reasonably possible. This means that the PR should be focused on a single change, bug fix, or feature addition. Multiple issues being address or unrelated fixes should be called out and addressed separately.
* There is no hard limit on the size of a PR, but less than 300 lines is a good target. However, this excludes documentation and tests.
* If an algorithm or expression is inefficient or could be improved, provide alternatives if possible. E.g. "This is too slow" is not helpful. "I think this can be done better by doing X instead of Y. What do you think?" is more helpful.

### Tests

* There should always be tests for customer-facing features and bug fixes. If there are no tests it should be explained in the PR description why there aren't.
* Tests should cover the "happy path" as well as edge cases and be as comprehensive as possible.

### Documentation

* Customer-facing features and bugfixes should include a note in `changelog.rst`.

### Security

* Check for any security vulnerabilities that might be introduced by this code.
* No private or sensitive information should be committed to the repository.
