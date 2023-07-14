function [milk,tea,filter,instant,sugar,total]=coffee_club(nmilk,ntea,nfilter,ninstant,nsugar)
%usage in number of [milk tea filter instant sugar] per day

milk_per_2litres=1.6; %pounds
milk_ml_per_cup=15; %1 teaspoon =5 ml

tea_cost = 3; %pounds per bag
nbags = 160;

filter_cost = 2.5; %pounds per bag
filter_grams = 227; %number of grams per bag of filter coffee

instant_cost = 6;
n_instant_grams = 300;

sugar_cost = 1;
n_sugar_grams = 500;

milk = milk_per_2litres*milk_ml_per_cup/2000;

tea = tea_cost/nbags;
n=7; filter = filter_cost/filter_grams*38 / n; %based on 38g (as measured) making n cups
instant = 1.7 * instant_cost / n_instant_grams; %measured 1.7g per spoon
sugar = 2*sugar_cost/n_sugar_grams;

total = ( nmilk*milk + ntea*tea + nfilter*filter + ninstant*instant + nsugar*sugar )*365/12;