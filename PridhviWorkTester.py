import random
import asyncio
from loguru import logger
from tqdm import tqdm

def generate_coeffs_answer_n(n=0):
    if n <= 0:
        n = random.randrange(100 + 1)
    sol = []
    for i in range(n):
        sol.append(random.randrange(-100, 100))
    # 		if i == 0:
    # 			sol[0] = 0
    # print(sol[-1])
    c = []
    a = []
    for i in range(n):
        coeffs = []
        for j in range(n):
            coeffs.append(random.randrange(-10, 10))
        # 			if i == 0 and j ==  0:
        # 				coeffs[0] = 0
        answer = 0
        for j in range(n):
            answer += coeffs[j] * sol[j]
        c.append(coeffs)
        a.append(answer)
    return c, a, n, sol


def write_trial(c, a, n, sol=None,dest=None):
    lines = []
    lines.append(str(n))
    for i in range(n):
        temp = []
        for j in c[i]:
            temp.append(str(j))
        lines.append(" ".join(temp))
    for i in range(n):
        lines.append(str(a[i]))
    if sol != None:
        lines.append("Solution")
        for l in sol:
            lines.append(str(l))
    if dest == None:
        for i in lines:
            print(i)
    else:
        with open(dest, "w") as f:
            for i in lines:
                f.write(str(i) + "\n")


# def write_trial_to_disk(n=4):
#     c, a, n, sol = generate_coeffs_answer_n(n)
#     write_trial(c, a, n, "data.dat")
#     return sol
fail_counter = 0
async def run_trial(n=4):
    global fail_counter
    c, a, n, sol = generate_coeffs_answer_n(n)
    # coeffs = str(n) + "\n"
    # coeffs = ["a.out <",coeffs]
    proc = await asyncio.subprocess.create_subprocess_exec("a.exe", stdin=asyncio.subprocess.PIPE, stdout=asyncio.subprocess.PIPE
    )
    nstr = str(n) + "\n"
    bytesn = bytes(nstr, encoding='utf-8')
    proc.stdin.write(bytesn)
    for s in c:
        for i in s:
            istr = str(i) + "\n"
            bytesi = bytes(istr, encoding='utf-8')
            proc.stdin.write(bytesi)
    for s in a:
        ss = str(s) + "\n"
        bs = bytes(ss, encoding='utf-8')
        proc.stdin.write(bs)
    # for i in range(1):
    #     proc.stdin.write(bytes("testing","utf-8"))
    # print(output)
    proc.stdin.drain()
    await proc.wait()
    for i in range(n):
        output = await proc.stdout.readline()
        output = output.decode("utf-8")
        output = output.strip()
        if(output!= str(sol[i])):
            fail_counter +=1
            # print(": " + str(a[i]))
            logger.error(f"Failed trial, wrote to disk.  {fail_counter} trials failed so far")
            write_trial(c,a,n,sol,"data.dat")
            return
    # logger.success("Finished Trial")
    

    # r = result.stdout.decode("utf-8")
    # r = r.strip().split()
    # if len(r) == len(sol):
    #     for answer in range(len(r)):
    #         if(r[answer]==str(sol[answer])):
    #             return True
    #         else:
    #             return c,a,n,sol
    # else:
    #     return c,a,n,sol

def run_trial_to_disk_if_fail(n=4):
    out = run_trial(n)
    if type(out) != bool:
        c,a,n,sol = out
        write_trial(c,a,n,"data.dat")
        print("Failed, dumped")
        print(" ".join(sol))

async def main():
    for i in tqdm(range(10000)):
        await run_trial(50)
loop = asyncio.get_event_loop()
loop.run_until_complete(main())
loop.close()
# main()