#!/bin/bash

for i in *.root.txt; do
    mv $i ${i/".root.txt"/".txt"}
done