# Shared SQL-identifier validation utilities

Internal helpers used by functions that accept a user-supplied
identifier (trait name, effect name, column name, model label) and then
interpolate it into SQL. Keeping the rules in one place guarantees every
caller enforces the same whitelist and rejects the same reserved
keywords.
