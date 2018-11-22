def find(arr, start_with=0):
    return [i+start_with for i in range(0, len(arr)) if arr[i]]

def greater(arr, num):
    return [x for x in arr if x>num]