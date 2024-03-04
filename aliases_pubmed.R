library(data.table)
library(stringr)
library(easyPubMed)

hgnc_complete_set <- fread('../data/hgnc_complete_set.txt', key = 'symbol')[
    !(ensembl_gene_id %in% names(table(ensembl_gene_id))[table(ensembl_gene_id) > 1]) & locus_group == 'protein-coding gene'
]

alias_table <- hgnc_complete_set[, .(symbol_alt = unique(c(str_split(alias_symbol, '\\|')[[1]], str_split(prev_symbol, '\\|')[[1]]))), by = symbol]
alias_table <- alias_table[symbol_alt != '']

# Remove elements of <symbol_alt> that also occur in <symbol>:
alias_table <- alias_table[!(symbol_alt %in% symbol)]

# Get number of PubMed search results for each alternative symbol:
alias_table[, pubmed_freq := as.numeric(get_pubmed_ids(paste0(symbol_alt, '[Text Word]'))$Count), by = symbol_alt]

# Adding combined searches, i.e. check for cases where the alias and the up-to-date symbol occur together:
alias_table[,
    pubmed_comb_freq := as.numeric(get_pubmed_ids(paste0('(', symbol_alt, '[Text Word]) AND (', symbol, '[Text Word])'))$Count),
    by = .(symbol, symbol_alt)
]

fwrite(alias_table, '../data/aliases_pubmed.csv')
