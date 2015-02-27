TESTOUT = toplevel_make_test.$$$$.output

tests:
	cd scripts/tests && make > $(TESTOUT) && diff $(TESTOUT) make.expected && rm -f $(TESTOUT)

.PHONY: tests
