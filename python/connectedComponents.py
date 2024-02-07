from PIL import Image
import numpy as np

def find_connected_components(image_path):
    # Load the bitmap image
    img = Image.open(image_path)
    img_array = np.array(img)

    def dfs(x, y, group_color):
        stack = [(x, y)]
        visited = set()
        component = []

        while stack:
            x, y = stack.pop()
            if (x, y) in visited or x < 0 or y < 0 or x >= img_array.shape[1] or y >= img_array.shape[0]:
                continue

            if img_array[y, x] == group_color:
                component.append((x, y))
                visited.add((x, y))
                stack.extend([(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)])

        return component

    connected_components = []
    group_color = 1

    for y in range(img_array.shape[0]):
        for x in range(img_array.shape[1]):
            if img_array[y, x] == 1:  # Assuming 1 represents connected pixels
                component = dfs(x, y, group_color)
                connected_components.append(component)
                group_color += 1

    return connected_components

image_path = "blob_test_2b.bmp"
components = find_connected_components(image_path)

for idx, component in enumerate(components):
    print(f"Connected Component {idx + 1}:")
    for x, y in component:
        print(f"Pixel at ({x}, {y})")

# You can also access the components as a list of lists:
# components[0] contains the coordinates of the first component, components[1] for the second, and so on.
