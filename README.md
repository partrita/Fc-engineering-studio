# Fc Engineering Studio

Fc_Hinge–CH2–CH3 영역을 파라미터화해 charge, glycosylation, FcγR binding hot-spot을 기반으로 설계하는 TUI(Text User Interface) 기반 도구입니다.

이 도구는 **EU Numbering** 체계를 기준으로 인간 IgG Fc 변이 서열을 신속하게 생성하고 검증할 수 있도록 도와줍니다.

## 주요 기능

- **Isotype & Allotype 선택**: IgG1, IgG2, IgG4의 다양한 Allotype(WT, Trastuzumab, Rituximab 등) 베이스 서열 지원
- **Preset Mutations (Common Mutations)**: 널리 알려진 주요 변이들을 체크박스로 간편하게 적용
  - LALA (Effector Silencing)
  - YTE/LS (Half-life extension)
  - Knob-into-Hole (Bispecific design)
  - S228P (IgG4 Hinge stabilization) 등
- **Custom Mutants 입력**: 사용자가 원하는 임의의 변이를 `S239D/I332E`와 같은 형식으로 직접 입력 가능
- **FASTA 생성 및 복사**: 생성된 서열을 즉시 FASTA 형식으로 확인하고 클립보드에 복사

## 설치 및 실행

이 프로젝트는 [uv](https://github.com/astral-sh/uv)를 사용하여 의존성을 관리합니다.

### 실행 방법
```bash
uv run fc-engineer
```

### 개발 환경 설정
```bash
# 의존성 설치
uv sync

# 테스트 실행
uv run pytest
```

## 사용법 및 단축키

1. 왼쪽 패널에서 **Isotype**과 **Allotype**을 선택합니다.
2. 적용하고자 하는 **Common Mutations**를 체크합니다.
3. 추가적인 변이가 있다면 **Custom Mutants** 창에 입력합니다 (슬래시`/` 또는 콤마`,`로 구분).
4. **Generate** 버튼을 누르거나 `Enter` 키를 입력하여 서열을 생성합니다.

### 단축키 (Shortcuts)
- `Enter`: 서열 생성 (Generate)
- `B`: 배치 모드 열기 (Mutation 화면에서)
- `Ctrl + Y`: 결과 서열 클립보드 복사 (Copy / Copy All)
- `Ctrl + G`: GenBank 형식으로 클립보드 복사
- `Ctrl + R`: CSV 리포트로 클립보드 복사
- `Ctrl + C`: 프로그램 종료 (Quit)

## 프로젝트 구조

```
src/fc_engineer/
├── __init__.py          # 패키지 메타데이터
├── __main__.py          # python -m fc_engineer 진입점
├── app.py               # TUI 화면 및 MutantApp
├── core.py              # EU Numbering 계산, 변이 파싱/적용, diff·배치·developability 로직
├── config.py            # YAML 데이터 로딩 및 검증
├── exporters.py         # GenBank / CSV 내보내기 포맷터
└── data/
    ├── sequences.yaml   # Isotype별 베이스 서열 데이터
    └── mutants.yaml     # 주요 변이 프리셋 데이터
tests/                   # 단위 테스트
pyproject.toml           # 프로젝트 설정 및 의존성 정의
```

## 데이터 커스터마이징

새로운 베이스 서열이나 변이 프리셋을 추가하려면 `src/fc_engineer/data/sequences.yaml` 또는 `src/fc_engineer/data/mutants.yaml` 파일을 수정하십시오. 프로그램 실행 시 해당 파일들을 자동으로 로드합니다.

## 로드맵 (Roadmap)

- [x] **IgG3 Isotype 추가**: IgG1, IgG2, IgG4에 이어 IgG3 베이스 서열 및 Allotype 지원 (확장 힌지 정렬 포함)
- [x] **변이 프리셋 확장**: 최신 문헌 기반의 developability 관련 변이(PTM hotspot 제거 등) 프리셋 추가
- [x] **서열 비교 뷰**: WT 대비 변이 서열의 위치별 diff 하이라이팅 (EU Numbering 기준)
- [x] **다양한 출력 포맷**: FASTA 외에 GenBank, CSV 리포트 등 내보내기 지원
- [x] **배치(Batch) 처리**: 여러 변이 조합을 일괄 생성하는 모드
- [x] **Developability 지표 요약**: glycosylation site(N297), charge 이슈, FcγR binding hotspot 주석 표시
