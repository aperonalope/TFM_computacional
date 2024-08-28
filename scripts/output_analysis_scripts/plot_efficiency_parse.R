setwd("/home/alvaro/Documentos/TFM/prueba_end/testeo")
library(dplyr)
library(ggplot2)
data <- read.csv("./simulaciones.txt", sep = " ", header =FALSE)[,2:5]
data[,3]=as.integer(data[,3])
colnames(data)= c('matches','read_name','Counts','Length')
args=commandArgs(trailingOnly = TRUE)
data <- read.csv(args[1], sep = " ", header =TRUE)
summed_data <- data %>%
  group_by(Length) %>% 
  summarise(total_counts = (sum(Counts))/10)

plot=ggplot(summed_data)+geom_point(aes(x=Length,y=total_counts))+geom_line(aes(x=Length,y=total_counts))+ylim(c(0,100))
plot



plot=ggplot(summed_data, aes(x = Length, y = total_counts)) +
    geom_line(color ="orange", size = 3) +
    ylim(80, 105) +
    scale_x_continuous(breaks = c(15, 50, 100, 135),expand = c(0.02, 0)) +
    theme_minimal(base_size = 15) +
    labs(x = "CircRNA length (nt)", y = "% of detected reads")+ theme(,axis.title = element_text(face = "bold"),
          axis.text = element_text(size = 16),
          panel.grid.major = element_line(color = "grey80", size = 0.5),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", color = "grey80"))
plot
write.csv(as.data.frame(summed_data),sep = " ",file = "./efficiency.csv")
ggsave("./plot.png", plot = plot, width = 7, height = 5, dpi = 300)


data_vicinal=read.csv("./resultados.txt",header = FALSE)
detected_vicinal=c()
lengths_vicinal=c()
for(i in seq(1,242)){
  
  if(i%%2==0){
    detected_vicinal=append(detected_vicinal,data_vicinal[i,1])
  }else{
    lengths_vicinal=append(lengths_vicinal,data_vicinal[i,1])
  }
}
detected_vicinal=detected_vicinal/10
tabla=cbind(summed_data[1:(length(row.names(summed_data))),],detected_vicinal)



library(ggplot2)

library(ggplot2)

plot = ggplot(tabla, aes(x = Length)) +
  geom_line(aes(y = total_counts, color = "Our pipeline"), size = 3) +
  geom_line(aes(y = detected_vicinal, color = "Vicinal 3.0"), size = 3) +
  scale_y_continuous(breaks = seq(0, 110, 10)) +
  scale_x_continuous(breaks = c(15, 50, 100, 135), expand = c(0.02, 0)) +
  theme_minimal(base_size = 15) +
  labs(x = "CircRNA length (nt)", y = "% of detected reads", color = "Software") +
  scale_color_manual(
    values = c("Our pipeline" = "orange", "Vicinal 3.0" = "blue"),
    name = "Software"  # Add the title "Software" to the legend
  ) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 16),
    panel.grid.major = element_line(color = "grey80", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "grey80")
  )


plot
ggsave("./plot_comparison.png", plot = plot, width = 10, height = 5, dpi = 300)
