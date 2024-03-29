# Shared server log in
# ssh root@10.19.139.174
# Biome944
# Save all data into data folder

# Import via manifest
  qiime tools import \
  --type "SampleData[SequencesWithQuality]" \
  --input-format SingleEndFastqManifestPhred33V2 \
  --input-path ./fish_ray_manifest.txt \
  --output-path ./demux.qza

# Creating a visualization of demultiplexed samples
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

# Denoise
 qiime dada2 denoise-single \
  --i-demultiplexed-seqs /data/SpinySamples/demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 235 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats.qza

  qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file /data/SpinySamples/Spiny_fish_metadata.txt

  qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

  screen -S denoise

   rep-seqs_spinyT3.qza

   # Generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs_spinyT3.qza \
  --o-alignment aligned-rep-seqs-SpinyT3.qza \
  --o-masked-alignment masked-aligned-rep-seqs-SpinyT3.qza \
  --o-tree unrooted-tree-SpinyT3.qza \
  --o-rooted-tree rooted-tree-SpinyT3.qza
