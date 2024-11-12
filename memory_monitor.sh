#!/bin/bash

LOG_FILE="memory_usage.log"
ACTIVE_TESTS_FILE="active_tests.json"

# Initialize the log file with headers
echo "Timestamp,MemoryUsage(MB),MemoryLimit(MB),ActiveTests" > "$LOG_FILE"

while true; do
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")

    # Read memory usage and limit
    MEMORY_USAGE=$(cat /sys/fs/cgroup/memory/memory.usage_in_bytes 2>/dev/null || echo 0)
    MEMORY_LIMIT=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || echo 0)
    MEMORY_USAGE_MB=$((MEMORY_USAGE / 1024 / 1024))
    MEMORY_LIMIT_MB=$((MEMORY_LIMIT / 1024 / 1024))

    # Read active tests
    if [ -f "$ACTIVE_TESTS_FILE" ]; then
        ACTIVE_TESTS=$(jq -r '.[]' "$ACTIVE_TESTS_FILE" | paste -sd "," -)
        if [ -z "$ACTIVE_TESTS" ]; then
            ACTIVE_TESTS="None"
        fi
    else
        ACTIVE_TESTS="No tests running"
    fi

    echo "$TIMESTAMP,$MEMORY_USAGE_MB,$MEMORY_LIMIT_MB,\"$ACTIVE_TESTS\"" >> "$LOG_FILE"
    sleep 10
done
