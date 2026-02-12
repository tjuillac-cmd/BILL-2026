library(ggplot2)

##df <- read.table("/home/totu/Documents/um_files/bill/diff-sniffles25-27.diff.sites_in_files", header=TRUE, sep="\t")
## df <- read.table("/home/totu/Documents/um_files/bill/diff-snp25-27.diff.sites_in_files", header=TRUE, sep="\t")
##df <- read.table("/home/totu/Documents/um_files/bill/diff-sniffles25-27cold.diff.sites_in_files", header=TRUE, sep="\t")
df <- read.table("/home/totu/Documents/um_files/bill/diff-snp25-27cold.diff.sites_in_files", header=TRUE, sep="\t")

# Distribution des variants par fichier
ggplot(df, aes(x=IN_FILE, fill=IN_FILE)) +
  geom_bar() +
  labs(title="Variants par fichier", x="Présent dans", y="Nombre de variants") +
  theme_minimal()

# Position des variants le long du génome
ggplot(df, aes(x=POS1, y=1, color=IN_FILE)) +
  geom_point(size=3) +
  labs(title="Distribution des variants sur la référence", x="Position (bp)") +
  theme_minimal()
