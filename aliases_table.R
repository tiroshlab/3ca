library(data.table)
library(magrittr)
library(matkot)

aliases_table <- fread('../data/aliases_pubmed.csv')
aliases_table <- aliases_table[
    pubmed_freq > 1000 &
        pubmed_comb_freq > 100 &
        !(symbol_alt %in% symbol) &
        str_count(symbol_alt, '[a-z]') <= 1 &
        str_count(symbol_alt, '.') > 1 & # Could put 2 here as some 2-letter aliases map poorly, but some are important like RB -> RB1
        !grepl('\\.', symbol_alt)
]
aliases_table <- aliases_table[!(symbol_alt %in% aliases_table[, .(l = length(unique(symbol))), by = symbol_alt][l > 1, symbol_alt])]
# Try capitalising and removing dashes. If this makes any of symbol_alt the same as symbol, remove these entries. Otherwise, if it makes any
# elements of symbol_alt the same, take the one with the highest pubmed_freq.
aliases_table[, symbol_alt_x := gsub('-', '', toupper(symbol_alt))]
aliases_table <- aliases_table[symbol_alt_x != symbol]
aliases_table <- aliases_table[,
    .(symbol_alt = rbindlist(lapply(
        names(table(symbol_alt_x)),
        function(symb) if(table(symbol_alt_x)[symb] == 1) {
            return(.SD[symbol_alt_x == symb, .(symbol_alt, pubmed_freq)])
        } else return(.SD[symbol_alt_x == symb, .(symbol_alt = symbol_alt[which.max(pubmed_freq)], pubmed_freq = max(pubmed_freq))])
    ))[order(-pubmed_freq), symbol_alt]),
    by = symbol
]

fwrite(aliases_table, '../data/alias_table.csv')
