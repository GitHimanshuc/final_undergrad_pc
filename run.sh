#!/bin/zsh
source ~/.zshrc
cd /home/himanshu/Desktop/final_year_now/new_2d_code_for_pc/
mpicc main.c -lm
mpiexec -np 4 a.out
python /home/himanshu/Desktop/final_year_now/new_2d_code_for_pc/just_ani2.py
