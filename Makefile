TESTOUT = toplevel_make_test.$$$$.output

tests:
	cd scripts/tests && make --no-print-directory > $(TESTOUT) && diff $(TESTOUT) make.expected && rm -f $(TESTOUT)

.PHONY: tests
