library(lubridate)
passofundo = c(63, 68, 35, 60, 55, 53, 13, 48, 26, 46, 93, 44, 85, 40, 34, 52, 
               16, 61, 7, 68, 71, 68, 17, 61, 58, 52, 64, 48, 33, 41, 38, 54, 
               85, 65, 41, 64, 27, 54, 60, 42, 72, 37, 51, 36, 79, 28, 65, 27, 
               74, 27, 74, 28, 39, 29, 84, 33, 56, 45, 69, 74, 76, 24, 63, 75, 
               49, 20, 39, 76, 32, 81, 31, 45, 37, 87, 42, 5, 44, 8, 48, 2, 59, 
               11, 57, 50, 44, 69, 43, 33, 52, 24, 52, 72, 45, 77, 35, 94, 28, 
               65, 31, 11, 44, 54, 71, 42, 99, 25, 94, 93, 92, 24, 82, 65, 89, 
               49, 81, 68, 75, 96, 77, 9)



data = paste(month(seq(ymd('2013-02-01'),ymd('2023-01-01'), by = '1 month'), 
                   label = TRUE), year(seq(ymd('2013-02-01'), ymd('2023-01-01'),
                                           by = '1 month')))

data_completa = seq(ymd('2018-02-01'),ymd('2023-01-01'), by = '1 month')

dados = data.frame('Data' = data,  'Capacidade Reservatório' = passofundo, 
                   'Data completa' = data_completa)

write.csv(dados, 'Capacidade_Reservatório.csv')
