import pytest
from fc_engineer.app import MutantApp, ResultScreen, BatchScreen
from textual.widgets import Input

@pytest.mark.asyncio
async def test_app_compose():
    app = MutantApp()
    async with app.run_test() as pilot:
        assert app.title == "Fc Engineering Studio"
        # If it reaches here, compose() didn't crash

@pytest.mark.asyncio
async def test_result_screen_renders_diff():
    app = MutantApp()
    async with app.run_test() as pilot:
        app.selected_isotype = "igg1"
        app.selected_allotype = "WT(P01857-1)"
        app.all_mutants = "L234A/L235A"
        await app.push_screen(ResultScreen())
        await pilot.pause()
        # FASTA 생성 및 diff 뷰 렌더링이 크래시 없이 완료되어야 함
        assert app.last_fasta.startswith(">IGG1_Wt(p01857-1)_L234A_L235A")
        mut_seq = app.last_fasta.split("\n")[1]
        assert mut_seq[116:118] == "AA"  # L234A/L235A 적용 확인 (idx 116=EU234)

@pytest.mark.asyncio
async def test_batch_screen_generation():
    app = MutantApp()
    async with app.run_test() as pilot:
        app.selected_isotype = "igg1"
        app.selected_allotype = "WT(P01857-1)"
        await app.push_screen(BatchScreen())
        await pilot.pause()
        screen = app.screen
        screen.query_one("#input-batch", Input).value = "L234A/L235A, N297A"
        screen.action_generate()
        await pilot.pause()
        assert len(app.batch_results) == 2
        headers = [h for h, _ in app.batch_results]
        assert headers[0] == "IGG1_Wt(p01857-1)_L234A_L235A"
        assert headers[1] == "IGG1_Wt(p01857-1)_N297A"

@pytest.mark.asyncio
async def test_batch_screen_back_wipes_state():
    app = MutantApp()
    async with app.run_test() as pilot:
        app.selected_isotype = "igg1"
        app.selected_allotype = "WT(P01857-1)"
        await app.push_screen(BatchScreen())
        await pilot.pause()
        app.screen.query_one("#input-batch", Input).value = "N297A"
        app.screen.action_generate()
        assert len(app.batch_results) == 1
        app.screen.action_back()
        await pilot.pause()
        # 뒤로 가기 시 배치 상태가 완전히 초기화되어야 함
        assert app.batch_results == []
