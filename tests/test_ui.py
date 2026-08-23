import pytest
from fc_engineer.app import MutantApp

@pytest.mark.asyncio
async def test_app_compose():
    app = MutantApp()
    async with app.run_test() as pilot:
        assert app.title == "Fc Engineering Studio"
        # If it reaches here, compose() didn't crash
