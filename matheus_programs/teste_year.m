close all;
clear all;


load ~/Dropbox/Mestrado/matheus_programs/PV_in_vars
cd 100/ 
load Tseries.mat                                   
load Tseries.mat                                   
junk_year = 2009:2014;
m_year = str2num(datestr(timeaxis,'yyyy'));  
month = str2num(datestr(timeaxis,'mm'));
m_date = [m_year';month'];
m_date = m_date';
num_of_years = length(junk_year);
        %%%% THE YEAR
        for j = 2%1:num_of_years;
                ind_year = find(m_year == junk_year(j))
                %%%% THE MONTH
                for k = 1:12;
                        ind_month = find(month(ind_year) == k)
                end;
        end;
cd ..
