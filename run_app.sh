#!/usr/bin/env bash
set -e
cd ~/Minor_Project_RNASeq
if [ -d "venv" ]; then
  source venv/bin/activate
fi
streamlit run app.py
