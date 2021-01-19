# mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
# model <- glm(formula = admit ~ gre * gpa + rank, data=mydata, family=binomial)
# summary(model)
# 
# 
# 
# interact_plot(model, pred = 'gpa' , modx = 'gre', plot.points = F, interval = T, modx.values = c(0,1), outcome.scale = 'link') 
#  
# values <- data.frame(gre = c(0,1), 
#                      gpa = c(3,3), 
#                      rank = rep(mean(mydata$rank, 2)))
# aa <- predict(model, values, type = 'link')
# aa
# 
# 
# 
# 
# values <- data.frame(gre = c(0,1), 
#                      gpa = c(4,4), 
#                      rank = rep(mean(mydata$rank, 2)))
# aa <- predict(model, values, type = 'link')
# aa
# 
# 
# interact_plot(model, pred = 'gpa' , modx = 'gre', plot.points = F, interval = T, modx.values = c(0,1), outcome.scale = 'response') 
# 
# 
# 
# 
# 
# # ok 
# values <- data.frame(gre = c(0,1), 
#                      gpa = c(1,1), 
#                      rank = rep(mean(mydata$rank, 2)))
# aa <- predict(model, values, type = 'link')
# aa <- exp(predict(model, values, type = 'link'))
# 2.362862e-05   / 2.330605e-05
# exp(0.018507-0.004762)
# 
# 
# 
# 
# # ok 
# values <- data.frame(gre = c(0,1), 
#                      gpa = c(3,3), 
#                      rank = rep(mean(mydata$rank, 2)))
# 
# aa <- exp(predict(model, values, type = 'link'))
# 0.03541974    / 0.03527052
# exp(0.018507-(3*0.004762))
# 
# 
# 
# # 
# values <- data.frame(gre = c(0,200), 
#                      gpa = c(3,3), 
#                      rank = rep(mean(mydata$rank, 2)))
# 
# aa <- exp(predict(model, values, type = 'link'))
# aa
# # 0.08205527    / 0.03527052
# 
# aa[2] / aa[1]
# 
# exp( (200*0.018507) - (200*3*0.004762))
# 
# 
# 
# exp( (500*0.018507) - (500*3*0.004762))
# exp( (250*0.018507) - (250*3*0.004762))
# 
# exp( (500*0.018507) - (500*3*0.004762) - (250*0.018507) + (250*3*0.004762)) # I think this is right answer... 
# 
# 
# 8.252366 / 2.872693
# 
# 
# 
# 
# tmp <- function(x) {
#   x / (1+x)
# }
# 
# 
# map_dbl(bb, tmp)
# 
# 
# 
# 
# 
