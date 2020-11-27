from loguru import logger
read_lines_file = []
logger.info("Attempting to read file.")
with open("data.dat") as f:
  read_lines_file = f.readlines()
  logger.success(f"Read {len(read_lines_file)} lines from file.")
length_matrix = int(read_lines_file[0].strip())
read_lines_file = read_lines_file[1:]
logger.success(f"Matrix size: {length_matrix}x{length_matrix}")
eqs = []
for j in range(length_matrix):
  l = read_lines_file[j]
  l = l.split()
  s = []
  for i in range(length_matrix):
    # logger.debug(l)
    s.append(l[i].strip() + "*" + str(chr(ord('a')+i)).strip())
  s = " + ".join(s)
  s += '='
  s += read_lines_file[length_matrix + j]
  s = s.strip()
  eqs.append(s)
print("; ".join(eqs))