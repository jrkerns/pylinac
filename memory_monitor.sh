#!/bin/bash

LOG_FILE="memory_usage.log"
echo "Timestamp,MemoryUsage(MB),MemoryLimit(MB)" > $LOG_FILE

while true; do
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
    MEMORY_USAGE=$(cat /sys/fs/cgroup/memory/memory.usage_in_bytes)
    MEMORY_LIMIT=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes)
    MEMORY_USAGE_MB=$((MEMORY_USAGE / 1024 / 1024))
    MEMORY_LIMIT_MB=$((MEMORY_LIMIT / 1024 / 1024))
    echo "$TIMESTAMP,$MEMORY_USAGE_MB,$MEMORY_LIMIT_MB" >> $LOG_FILE
    sleep 10
done
