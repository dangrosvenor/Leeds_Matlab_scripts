function [days_required_for_mean,time_mean_str] = days_of_month(month_no)

            switch month_no
                case 1
                    days_required_for_mean = [1:31]; time_mean_str = 'Jan';
                case 2
                    days_required_for_mean = [32:60]; time_mean_str = 'Feb';
                case 3
                    days_required_for_mean = [61:91]; time_mean_str = 'Mar';
                case 4
                    days_required_for_mean = [92:121]; time_mean_str = 'Apr';
                case 5
                    days_required_for_mean = [122:152]; time_mean_str = 'May';
                case 6
                    days_required_for_mean = [153:182]; time_mean_str = 'Jun';
                case 7
                    days_required_for_mean = [183:213]; time_mean_str = 'Jul';
                case 8
                    days_required_for_mean = [214:244]; time_mean_str = 'Aug';
                case 9
                    days_required_for_mean = [245:274]; time_mean_str = 'Sep';
                case 10
                    days_required_for_mean = [275:305]; time_mean_str = 'Oct';
                case 11
                    days_required_for_mean = [306:335]; time_mean_str = 'Nov';
                case 12
                    days_required_for_mean = [336:366]; time_mean_str = 'Dec';
            end
            