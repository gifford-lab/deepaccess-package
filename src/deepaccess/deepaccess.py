#!/usr/bin/env python
# -*- coding: utf-8 -*-
import deepaccess.train as train
import deepaccess.interpret as interpret


def main():
    command_list = [train, interpret]
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=[command_list])
    parser.add_argument("command_args",
                        help="command arguments",
                        nargs=argparse.REMAINDER
                        )
    opts = parser.parse_args()
    command = opts.command
    command_args = opts.command_args
    args = parser.parse_args(command_args)
    command.main(args)
    
if __name__ == "__main__":
    main()
