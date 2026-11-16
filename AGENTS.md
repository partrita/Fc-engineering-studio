# Fc Engineering Studio - Project Context

## Project Overview
Fc Engineering Studio는 인간 IgG 및 이중특이성 항체(Monospecific IgG, CrossMab, IgG-scFv, scFv-Fc, DVD-Ig, BiTE 등)의 **도메인 토폴로지**와 **EU Numbering** 기반 변이 서열을 시각적으로 설계하고 분석하는 Python 기반 Flet 데스크톱 애플리케이션입니다.

### 핵심 기술 스택
- **Language**: Python >= 3.14
- **GUI Framework**: [Flet](https://flet.dev/) (`flet[desktop]>=0.86.0`)
- **Image Export**: [Pillow](https://python-pillow.org/) (`pillow>=11.0.0`, 투명 RGBA `.png` 모식도 렌더링)
- **Dependency Management**: [uv](https://github.com/astral-sh/uv)
- **Data Persistence**: YAML (`pyyaml`), GenBank-style `.ab` 프로젝트 파일, multi-FASTA, CSV
- **Clipboard Management**: `pyperclip` (30초 자동 삭제 보안 타이머 + 종료 시 정리)
- **Testing**: `pytest`, `pytest-asyncio`

### 아키텍처 (src-layout 패키지)
- `src/fc_engineer/app.py`: Flet 데스크톱 UI(`build_studio_ui`, `main`)와 상태 컨트롤러(`AntibodyStudioController`)를 포함합니다.
- `src/fc_engineer/designer.py`: 모듈형 도메인/체인 모델(`DomainBlock`, `PolypeptideChain`, `InterChainPairing`, `AntibodyProject`), 7종 포맷 템플릿(`build_format_template`), 원자적 교차 도메인 변이 라우팅(`route_mutations_to_chain`), 2D 모식도 씬 모델(`build_schematic_scene`)을 제공합니다.
- `src/fc_engineer/core.py`: EU Numbering 인덱스 계산(`get_residue_index`, `get_eu_position`), 상수 도메인 슬라이싱(`slice_constant_domain`, `get_domain_residue_index`, `get_domain_eu_position`), 변이 파싱/적용(`parse_mutation_list`, `apply_mutations`, `apply_domain_mutations`), diff 및 developability 분석 순수 로직을 포함합니다.
- `src/fc_engineer/config.py`: YAML 데이터 로딩(`load_yaml_data`, `load_domain_templates`, `_load_sequences_bundle`)과 보안 검증(`NoAliasSafeLoader`, 1MB 파일 크기 제한 `MAX_FILE_SIZE`)을 담당합니다.
- `src/fc_engineer/exporters.py`: `.ab` 저장/파싱(`serialize_ab_project`, `parse_ab_project`, `save_ab_file`, `load_ab_file`), multi-FASTA, GenBank, CSV 및 Pillow 기반 투명 PNG 내보내기(`export_schematic_png`, `save_schematic_png`)를 제공합니다.
- `src/fc_engineer/__main__.py`: `python -m fc_engineer` 실행을 지원합니다.
- `src/fc_engineer/data/sequences.yaml`: IgG1–4 Isotype/Allotype 상수 영역과 기본 `VH`/`VL`/`CL`/`G4S` 템플릿 서열을 관리합니다.
- `src/fc_engineer/data/mutants.yaml`: 읽기 전용 기본 Fc 변이 프리셋을 정의합니다.
- `tests/`: 도메인 EU 번호 체계, 교차 도메인 라우팅, `.ab`/`.png`/FASTA 내보내기 및 보안 검증, Flet 컨트롤러 및 헤드리스 UI 테스트를 포함합니다.

## Building and Running

### 실행 명령어
```bash
# 애플리케이션 실행
uv run fc-engineer

# 또는 모듈로 직접 실행
uv run python -m fc_engineer
```

### 개발 및 테스트
```bash
# 의존성 동기화
uv sync

# 테스트 실행
uv run pytest
```

## Development Conventions

### Coding Style
- **Python Typing**: `typing` 모듈을 사용하여 적극적으로 타입을 명시합니다 (`Dict`, `List`, `Optional`, `Tuple` 등).
- **Docstrings**: 클래스 및 주요 함수에 대한 설명을 포함하여 코드 가독성을 유지합니다.
- **Entry Points**: `[project.scripts]`를 통해 터미널 명령어를 제공하며, `main()` 함수를 두어 실행을 래핑합니다.

### Testing Practices
- `pytest` 및 `pytest-asyncio`를 사용하여 순수 도메인/EU 로직, 파일 직렬화·보안 검증, 컨트롤러/헤드리스 Flet UI 이벤트를 검증합니다.
- 특히 상수 도메인 EU 슬라이싱, 교차 도메인 원자적 롤백, `.ab` 파싱 보안 경계, 클립보드 타이머 경쟁 상태 테스트를 중시합니다.

### Data Management
- 패키지 기본 서열과 프리셋은 `src/fc_engineer/data/*.yaml`에 선언적으로 관리하며 읽기 전용으로 취급합니다 (`pyproject.toml`의 `[tool.setuptools.package-data]` 포함).
- 사용자가 런타임에 추가하는 커스텀 변이는 `mutants.yaml`을 수정하지 않고 현재 `.ab` 프로젝트 파일의 `file_custom_mutations`에만 저장합니다.

## User Interface Guidelines
- **Toolbar & DnD Interactions**: 상단 툴바 액션(프로젝트 이름, 포맷 선택, `Add Chain`, `Load .ab`, `Save .ab`, `Export FASTA`, `Export PNG`, `Copy FASTA`)과 좌측 팔레트 드래그 앤 드롭(`ft.Draggable` / `ft.DragTarget`), 2D 캔버스 드래그(`ft.GestureDetector`)를 통해 직관적으로 조립합니다.
- **Feedback & Lifecycle**: 오류 및 작업 완료 알림은 `page.show_dialog(ft.SnackBar(...))`로 표시하고, 창 닫기 이벤트(`ft.WindowEventType.CLOSE`) 및 `atexit` 시 `controller.shutdown()`을 호출하여 30초 클립보드 타이머 및 OS 클립보드 잔여 데이터를 정리합니다.

