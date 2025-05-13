
with open('./data/etopo1to5.bin', 'rb') as f:
    bytes = f.read(2)
    res = int.from_bytes(bytes, 'little', signed=True)
    print(res)
