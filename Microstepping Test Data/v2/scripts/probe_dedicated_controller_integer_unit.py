#!/usr/bin/env python3
'''Safely identify the physical unit of one integer MR position count.'''

import time

import serial


PORT = 'COM5'
BAUD = 115200
AXIS = 'X'
COMMAND_TIMEOUT_S = 5.0
STATUS_TIMEOUT_S = 30.0


def command(link, text, expected_prefix=None):
    expected = expected_prefix or text
    link.reset_input_buffer()
    print('SEND:', text)
    link.write((text + '\r').encode('ascii'))
    link.flush()
    deadline = time.time() + COMMAND_TIMEOUT_S
    received = []
    while time.time() < deadline:
        raw = link.readline()
        if not raw:
            continue
        line = raw.decode('ascii', errors='replace').strip()
        if not line:
            continue
        received.append(line)
        print('RECEIVED:', line)
        if line == 'COMM_ERR':
            raise RuntimeError('Controller rejected: ' + text)
        if line.upper().startswith(expected.upper()):
            return line
    raise RuntimeError(
        'No expected response for {0}; received {1}'.format(text, received)
    )


def wait_ready(link):
    deadline = time.time() + STATUS_TIMEOUT_S
    while time.time() < deadline:
        response = command(link, 'DS X', 'DS X ')
        status = int(response.split()[-1])
        if status == 1:
            return
        if status != 2:
            raise RuntimeError('Unexpected motor status: ' + response)
        time.sleep(0.05)
    raise RuntimeError('Motor did not become ready')


def main():
    link = serial.Serial(
        port=PORT, baudrate=BAUD, bytesize=8, parity='N', stopbits=1,
        timeout=0.2, write_timeout=2.0,
    )
    try:
        command(link, 'DM', 'DM ')
        command(link, 'SM X 4 200 0')
        command(link, 'SC X 400 400')
        command(link, 'SS X 471 471 628 628 1')
        command(link, 'ME X')
        wait_ready(link)
        command(link, 'SP X 0')

        print()
        input('Start IDS recording, then press Enter to command MR X 1. ')
        command(link, 'MR X 1', 'MR X ')
        wait_ready(link)
        command(link, 'DP X', 'DP X ')
        time.sleep(3.0)

        command(link, 'MR X -1', 'MR X ')
        wait_ready(link)
        final_position = command(link, 'DP X', 'DP X ')
        time.sleep(3.0)
        print()
        print('Returned position:', final_position)
        input('Stop IDS recording, then press Enter to disable the motor. ')
    finally:
        try:
            command(link, 'MO X', 'MO X')
        except Exception as exc:
            print('WARNING: MO X failed:', exc)
        link.close()


if __name__ == '__main__':
    main()
