import pandas as pd


class RegionIterator:
    """Iterate over the regions in the database"""

    def __init__(self, regions: pd.DataFrame) -> None:
        self.regions = regions.iterrows()

    def __aiter__(self) -> "RegionIterator":
        return self

    async def __anext__(self) -> tuple[str, pd.Series]:
        try:
            region = next(self.regions)
        except StopIteration:
            raise StopAsyncIteration
        return region
