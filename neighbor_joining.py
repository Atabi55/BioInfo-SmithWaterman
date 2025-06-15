inputs = input("Enter inputs: ")
inputs = inputs.split()
n = len(inputs)

matrix = []
distance = [0] * n
Q = [[0]*n for _ in range(n)]

for x in range(n):
    row = input("Enter values: ")
    values = list(map(int,row.split()))
    matrix.append(values)

print("\n My Matrix:")
for row in matrix:
    print(row)

for x in range(n):
    total = 0
    for y in range(n):
        if (y != x):
            distance[x] += matrix[x][y]

print("\nOur Distance:")
for x in distance:
    print(x)

print("\nQ matrix:")
for x in range(n):
    for y in range(n):
        if (y != x):
             Q[x][y] = (n-2)*matrix[x][y] - distance[x] - distance[y]

for row in Q:
    print(row)