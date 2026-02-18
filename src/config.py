import os


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB_NAME = "immune_cells.db"
DB_PATH = os.path.join(ROOT_DIR, DB_NAME)
CSV_FILE = os.path.join(ROOT_DIR, "cell-count.csv")


CELL_TYPES = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
