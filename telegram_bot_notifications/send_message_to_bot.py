import telegram_send
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-m', '--message', required=True,
    help="Message to send")
args = parser.parse_args()

telegram_send.send(messages=[f"{args.message}"])
