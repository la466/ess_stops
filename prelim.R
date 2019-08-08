library(ggplot2)

# just stop codon density either side of the splice site
file = read.csv("results/temp/stop_codon_position_densities.csv", head = T)
cairo_pdf("results/plots/all_stop_density.pdf")
plot(file$position, file$density, xlab = "Position relative to splice site", ylab = "Stop codon density")
abline(v = 0, lty = 2)
dev.off()

# stop density in each frame each side of the splice site
file = read.csv("results/temp/stop_codon_position_frame_densities.csv", head = T)
file$frame = as.factor(file$frame)
plot = ggplot(data = file, aes(x = file$position, y = file$density, col = file$frame)) + 
  geom_line() + 
  geom_vline(xintercept = 0, lty = 2)
ggsave(plot = plot, filename = "results/plots/stop_density_per_frame.pdf", width= 10, height = 5)

# stop density in each frame each side of the splice site for 3n introns
file = read.csv("results/temp/stop_codon_position_frame_3n_densities.csv", head = T)
file$frame = as.factor(file$frame)
plot = ggplot(data = file, aes(x = file$position, y = file$density, col = file$frame)) + 
  geom_line() + 
  geom_vline(xintercept = 0, lty = 2)
ggsave(plot = plot, filename = "results/plots/stop_density_per_frame_3n_introns.pdf", width= 10, height = 5)


# stop density in each frame at the first stop codon hit
file = read.csv("results/temp/stop_codon_position_frame_densities_first_hit.csv", head = T)
file$frame = as.factor(file$frame)
plot = ggplot(data = file, aes(x = file$position, y = file$density, col = file$frame)) + 
  geom_line() + 
  geom_vline(xintercept = 0, lty = 2)
ggsave(plot = plot, filename = "results/plots/stop_density_per_frame_3n_introns_first_hit.pdf", width= 10, height = 5)

# cassete exons
file = read.csv("results/temp/cassette_frame_densities.csv", head = T)
file$frame = as.factor(file$frame)
plot = ggplot(data = file, aes(x = file$position, y = file$density, col = file$frame)) + 
  geom_line() + 
  geom_vline(xintercept = 0, lty = 2)
plot
ggsave(plot = plot, filename = "results/plots/casstte.pdf", width= 10, height = 5)

# cassete exon decoys
file = read.csv("results/temp/cassette_decoys_frame_densities.csv", head = T)
file$frame = as.factor(file$frame)
plot = ggplot(data = file, aes(x = file$position, y = file$density, col = file$frame)) + 
  geom_line() + 
  geom_vline(xintercept = 0, lty = 2)
plot
ggsave(plot = plot, filename = "results/plots/casstte_decoys.pdf", width= 10, height = 5)

# cassete exon decoys 3n
file = read.csv("results/temp/cassette_decoys_3n_frame_densities.csv", head = T)
file$frame = as.factor(file$frame)
plot = ggplot(data = file, aes(x = file$position, y = file$density, col = file$frame)) + 
  geom_line() + 
  geom_vline(xintercept = 0, lty = 2)
plot
ggsave(plot = plot, filename = "results/plots/casstte_decoys_3n.pdf", width= 10, height = 5)

# constitutive exons
file = read.csv("results/temp/constitutive_frame_densities.csv", head = T)
file$frame = as.factor(file$frame)
plot = ggplot(data = file, aes(x = file$position, y = file$density, col = file$frame)) + 
  geom_line() + 
  geom_vline(xintercept = 0, lty = 2)
plot
ggsave(plot = plot, filename = "results/plots/constitutive.pdf", width= 10, height = 5)


# constitutive exons 3n
file = read.csv("results/temp/constitutive_3n_frame_densities.csv", head = T)
file$frame = as.factor(file$frame)
plot = ggplot(data = file, aes(x = file$position, y = file$density, col = file$frame)) + 
  geom_line() + 
  geom_vline(xintercept = 0, lty = 2)
plot
ggsave(plot = plot, filename = "results/plots/constitutive_3n.pdf", width= 10, height = 5)

# constitutive decoys
file = read.csv("results/temp/constitutive_decoy_frame_densities.csv", head = T)
file$frame = as.factor(file$frame)
plot = ggplot(data = file, aes(x = file$position, y = file$density, col = file$frame)) + 
  geom_line() + 
  geom_vline(xintercept = 0, lty = 2)
plot
ggsave(plot = plot, filename = "results/plots/constitutive_decoy.pdf", width= 10, height = 5)


# constitutive decoys 3n
file = read.csv("results/temp/constitutive_decoy_3n_frame_densities.csv", head = T)
file$frame = as.factor(file$frame)
plot = ggplot(data = file, aes(x = file$position, y = file$density, col = file$frame)) + 
  geom_line() + 
  geom_vline(xintercept = 0, lty = 2)
plot
ggsave(plot = plot, filename = "results/plots/constitutive_decoy_3n.pdf", width= 10, height = 5)



