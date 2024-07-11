def remove_overlapping_entities(train_data):
    corrected_data = []
    for text, annotations in train_data:
        entities = annotations['entities']
        # Sort entities by start position
        entities = sorted(entities, key=lambda x: x[0])

        # Initialize variables for the corrected entities
        corrected_entities = []
        prev_start, prev_end = None, None

        for start, end, label in entities:
            if prev_start is None:
                prev_start, prev_end = start, end
                corrected_entities.append((start, end, label))
            else:
                # Check for overlap
                if start < prev_end:
                    # Overlapping entity found, resolve conflict (e.g., skip or adjust)
                    if end > prev_end:
                        # Choose the longer entity or the one with higher priority
                        # For simplicity, we'll skip overlapping entities in this example
                        print(f"Skipping overlapping entity: ({start}, {end}, {label})")
                    continue
                else:
                    # No overlap, add to corrected entities
                    corrected_entities.append((start, end, label))
                    prev_start, prev_end = start, end

        corrected_data.append((text, {'entities': corrected_entities}))

    return corrected_data


# Example usage
TRAIN_DATA = [
    ('text1', {'entities': [(106, 107, 'ACTIVITY'), (110, 120, 'RECEPTOR')]}),
    ('text2', {'entities': [(100, 105, 'ACTIVITY'), (102, 112, 'RECEPTOR')]})
    # Add more data as needed
]

# Remove overlapping entities
corrected_data = remove_overlapping_entities(TRAIN_DATA)

# Print the corrected data (for verification)
for text, annotations in corrected_data:
    print(text)
    print(annotations)

# Now you can proceed with your SpaCy training using the corrected_data
