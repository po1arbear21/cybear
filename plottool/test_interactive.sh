#!/bin/bash
# Test script for interactive mode - demonstrates field checking

echo "Testing interactive mode with field checking..."
echo ""

# Test input:
# 1 - Select first job (efield_test)
# 1 - Select first file (usually device.fbs or similar)
# 3 - Check available fields
# a - Auto mode
# . - Current directory

echo -e "1\n1\n3\na\n.\ny" | python3 interactive.py 2>/dev/null

echo ""
echo "Test complete! Check the output above."