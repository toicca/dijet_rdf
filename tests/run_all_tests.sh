#!/bin/bash

for test_script in test/*.sh; do
  if [ "$test_script" != "test/run_all_tests.sh" ]; then
    bash "$test_script"
    RET=$?
    if [ $RET -ne 0 ]; then
      echo "Test $test_script failed with exit code $RET"
      exit $RET
    fi
  fi
done

echo "All tests passed successfully."
exit 0
