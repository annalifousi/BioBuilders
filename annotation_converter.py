import pandas as pd
import json
import os


with open('annotations.json', 'r') as f:
    data = json.load(f)

entity_name = "ACTIVITY"

train_data = data['annotations']
train_data = [tuple(i) for i in train_data]
for i in train_data:
    if i[1]['entities'] == []:
        i[1]['entities'] = (0, 0, entity_name)
    else:
        i[1]['entities'][0] = tuple(i[1]['entities'][0])

print(train_data)