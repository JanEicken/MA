#!/bin/bash
echo "Moving Log"
mv backend.log backend-$(date "+%Y%m%d%H%M%S").log
echo "Starting Worker"
python backend/backend-server.py 2>&1 >> backend.log &
PID=$!
echo Worker [$PID]
python server.py
echo "Stopping Worker"
kill $PID

