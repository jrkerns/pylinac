#!/bin/bash

# Configuration
LOG_FILE="memory_usage.log"
ACTIVE_TESTS_FILE="active_tests.json"
# Dynamically set TEMP_DIR based on TMPDIR environment variable; default to /tmp if TMPDIR is not set
TEMP_DIR="${TMPDIR:-/tmp}"
POLL_INTERVAL=10  # Polling interval in seconds

# Initialize the log file with headers
echo "Timestamp,MemoryUsage(MB),MemoryLimit(MB),ActiveTests,TempDirSize(MB)" > "$LOG_FILE"

while true; do
    # Capture the current timestamp
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")

    # Read memory usage and limit
    if [ -f /sys/fs/cgroup/memory/memory.usage_in_bytes ] && [ -f /sys/fs/cgroup/memory/memory.limit_in_bytes ]; then
        MEMORY_USAGE=$(cat /sys/fs/cgroup/memory/memory.usage_in_bytes 2>/dev/null || echo 0)
        MEMORY_LIMIT=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes 2>/dev/null || echo 0)
    else
        # Fallback for systems without cgroup memory files
        MEMORY_USAGE=$(free -b | awk '/Mem:/ {print $3}' || echo 0)
        MEMORY_LIMIT=$(free -b | awk '/Mem:/ {print $2}' || echo 0)
    fi

    MEMORY_USAGE_MB=$((MEMORY_USAGE / 1024 / 1024))
    MEMORY_LIMIT_MB=$((MEMORY_LIMIT / 1024 / 1024))

    # Read active tests
    if [ -f "$ACTIVE_TESTS_FILE" ]; then
        # Use jq to parse the JSON array and concatenate test names
        ACTIVE_TESTS=$(jq -r '.[]' "$ACTIVE_TESTS_FILE" | paste -sd "," -)
        # Handle empty active tests
        if [ -z "$ACTIVE_TESTS" ]; then
            ACTIVE_TESTS="None"
        fi
    else
        ACTIVE_TESTS="No tests running"
    fi

    # Determine the temporary directory to monitor
    # Prefer TMPDIR if set; else default to /tmp
    if [ -n "$TMPDIR" ]; then
        CURRENT_TEMP_DIR="$TMPDIR"
    else
        CURRENT_TEMP_DIR="/tmp"
    fi

    # Calculate the size of the temporary directory in MB
    if [ -d "$CURRENT_TEMP_DIR" ]; then
        # Use du to calculate the size. Suppress errors for directories with restricted permissions.
        TEMP_DIR_SIZE=$(du -sm "$CURRENT_TEMP_DIR" 2>/dev/null | awk '{print $1}')
        # Handle cases where du fails
        if [ -z "$TEMP_DIR_SIZE" ]; then
            TEMP_DIR_SIZE="Unknown"
        fi
    else
        TEMP_DIR_SIZE="Directory not found"
    fi

    # Log the data
    echo "$TIMESTAMP,$MEMORY_USAGE_MB,$MEMORY_LIMIT_MB,\"$ACTIVE_TESTS\",$TEMP_DIR_SIZE" >> "$LOG_FILE"

    # Wait for the next poll
    sleep "$POLL_INTERVAL"
done
