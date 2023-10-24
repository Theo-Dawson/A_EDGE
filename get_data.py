import yfinance as yf

msft = yf.Ticker("MSFT")

# get all stock info
print(msft.info)

# get historical market data
hist = msft.history(period="1mo")

# show meta information about the history (requires history() to be called first)
print(msft.history_metadata)

# show actions (dividends, splits, capital gains)
print(msft.actions)
print(msft.dividends)
print(msft.splits)
print(msft.capital_gains)  # only for mutual funds & etfs

# show share count
msft.get_shares_full(start="2022-01-01", end=None)

# show financials:
# - income statement
print(msft.income_stmt)
print(msft.quarterly_income_stmt)
# - balance sheet
print(msft.balance_sheet)
print(msft.quarterly_balance_sheet)
# - cash flow statement
print(msft.cashflow)
print(msft.quarterly_cashflow)
# see `Ticker.get_income_stmt()` for more options

# show holders
print(msft.major_holders)
print(msft.institutional_holders)
print(msft.mutualfund_holders)

# Show future and historic earnings dates, returns at most next 4 quarters and last 8 quarters by default. 
# Note: If more are needed use msft.get_earnings_dates(limit=XX) with increased limit argument.
print(msft.earnings_dates)

# show ISIN code - *experimental*
# ISIN = International Securities Identification Number
print(msft.isin)

# show options expirations
print(msft.options)

# show news
print(msft.news)

# get option chain for specific expiration
opt = msft.option_chain('YYYY-MM-DD')
# data available via: opt.calls, opt.puts


## PROXY SERVER FOR DOWNLOADING DATA EXAMPLE
# msft = yf.Ticker("MSFT")

# msft.history(..., proxy="PROXY_SERVER")
# msft.get_actions(proxy="PROXY_SERVER")
# msft.get_dividends(proxy="PROXY_SERVER")
# msft.get_splits(proxy="PROXY_SERVER")
# msft.get_capital_gains(proxy="PROXY_SERVER")
# msft.get_balance_sheet(proxy="PROXY_SERVER")
# msft.get_cashflow(proxy="PROXY_SERVER")
# msft.option_chain(..., proxy="PROXY_SERVER")

## MULTIPLE TICKET SYMBOLS
# tickers = yf.Tickers('msft aapl goog')

# # access each ticker using (example)
# tickers.tickers['MSFT'].info
# tickers.tickers['AAPL'].history(period="1mo")
# tickers.tickers['GOOG'].actions

##DOWNLOAD PRICE HISTORY INTO ONE TABLE
# data = yf.download("SPY AAPL", period="1mo")