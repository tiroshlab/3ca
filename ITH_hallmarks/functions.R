make_alias_table <- function(hgnc_complete_set) {
    out <- hgnc_complete_set[,
        rbind(
            data.table(symbol_alt = str_split(alias_symbol, '\\|')[[1]], type = 'alias'),
            data.table(symbol_alt = str_split(prev_symbol, '\\|')[[1]], type = 'prev')
        ),
        by = symbol
    ][symbol_alt == '', symbol_alt := NA]
    setkey(out, symbol_alt)
    return(out)
}

update_symbols_fast <- function(symbols, alias_table) {
    if(length(symbols) != length(unique(symbols))) warning('Not all symbols are unique!')
    setkey(alias_table, symbol_alt)
    symbols[grep('C[0-9]+ORF[0-9]+', symbols)] <- gsub('ORF', 'orf', symbols[grep('C[0-9]+ORF[0-9]+', symbols)])
    cond <- !(symbols %in% alias_table$symbol)
    # The following addition to <cond> ensures two different symbols don't map to the same "prev" symbol...
    cond <- cond & !(symbols %in% alias_table[type == 'prev'][symbols[cond]][symbol %in% names(table(symbol))[table(symbol) > 1], symbol_alt])
    # ...And the following makes sure we have all unique prev symbols:
    cond <- cond & !(symbols %in% alias_table[type == 'prev'][symbols[cond], names(table(symbol_alt))[table(symbol_alt) > 1]])
    if(any(cond)) {
        # alias_table[type == 'prev'][symbols[cond]]
        symbols[cond] <- alias_table[type == 'prev'][symbols[cond]][,
            .(newsymb = { # Match symbols that have at least one "prev" symbol that isn't already in symbols:
                cands <- symbol[!is.na(symbol) & !(symbol %in% symbols)]
                ifelse(length(cands) == 1, cands, symbol_alt) # Convert only those symbols that have exactly one such "prev" symbol
            }),
            by = symbol_alt
        ]$newsymb
    }
    cond <- !(symbols %in% alias_table$symbol)
    # As above, the following additions to <cond> ensure two different symbols don't map to the same "alias", and we have all unique aliases:
    cond <- cond & !(symbols %in% alias_table[type == 'alias'][symbols[cond]][symbol %in% names(table(symbol))[table(symbol) > 1], symbol_alt])
    cond <- cond & !(symbols %in% alias_table[type == 'alias'][symbols[cond], names(table(symbol_alt))[table(symbol_alt) > 1]])
    if(any(cond)) {
        symbols[cond] <- alias_table[type == 'alias'][symbols[cond]][,
            .(newsymb = { # Match symbols that have at least one alias that isn't already in symbols:
                cands <- symbol[!is.na(symbol) & !(symbol %in% symbols)]
                ifelse(length(cands) == 1, cands, symbol_alt) # Convert only those symbols that have exactly one such alias
            }),
            by = symbol_alt
        ]$newsymb
    }
    return(symbols)
}

dirr <- function(x) dir(x, recursive = TRUE)
