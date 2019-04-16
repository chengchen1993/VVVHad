#!/bin/bash
crab submit -c crab3_analysisQCD100-200.py
crab submit -c crab3_analysisQCD200-300.py
crab submit -c crab3_analysisQCD300-500.py
crab submit -c crab3_analysisQCD500-700.py
crab submit -c crab3_analysisQCD700-1000.py
crab submit -c crab3_analysisQCD1000-1500.py
crab submit -c crab3_analysisQCD1500-2000.py
crab submit -c crab3_analysisQCD2000-Inf.py


