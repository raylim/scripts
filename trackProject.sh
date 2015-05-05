#!/usr/bin/env bash

find . -ls -not -name '.*' > fls.txt 
git add fls.txt

git add sample_sets.txt
git add samples.txt
stat -t *.bed &> /dev/null && git add *.bed
stat -t alltables/*.txt &> /dev/null && git add alltables/*.txt
stat -t metrics/* &> /dev/null && git add metrics/*
stat -t version/* &> /dev/null && git add version/*

git commit -m "$(date) : project update"
git push
