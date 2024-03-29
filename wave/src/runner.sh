#!/bin/bash

for Nxt in  1000 # 500 1000 5000 10000 15000
	do 
		echo "#--------------| Nxt = $Nxt |---------------"
		(echo "bcNonRUpwindO2 ic1 GRk2 LeapFrog"
		echo "0 0  0"						# for bc :
		echo "0 5 6"						# for ic : mean sigma amplitude 
		echo "2 -3 1 0"						# for pot: e.g. (PT): V0, a    , x0
											#               (GR): l , sigma, rs
		echo "150 -80 80"      				# tmax, x*min, x*max
		echo "$Nxt" 
		echo "0.8"      )| python3 ./main.py			# Nt , Nx 
		echo "#-------------------------------------------"
	done
#----------------------------------------------------------------------
#                         .////.
#                         /// 6|
#                         //  _|
#                        _/_,-'
#                   _.-/'/   \   ,/;,
#                ,-' /'  \_   \ / _/
#                `\ /     _/\  ` /
#                  |     /,  `\_/
#                  |      \'
#      /\_        /`      /\
#    /' /_``--.__/\  `,. /  \
#   |_/`  `-._     `\/  `\   `.
#             `-.__/'     `\   |
#                           `\  \
#                             `\ \
#                               \_\__
#                                \___)
#---------------------------------------------------------------------