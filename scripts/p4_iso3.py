"""ISO 3166-1 alpha-3 codes to the country names used in p4_build_groundtruth.py (only countries in our region tables).
Used for structured names such as RVA/Cow-tc/NGA/NGRBg8/1998/G8P[1], where the 3-letter code is the independent place evidence."""
ISO3 = {
    # Africa
    "KEN": "kenya", "UGA": "uganda", "NGA": "nigeria", "GHA": "ghana", "ZAF": "south africa", "EGY": "egypt", "ETH": "ethiopia", "TZA": "tanzania", "CMR": "cameroon",
    "SEN": "senegal", "MLI": "mali", "ZMB": "zambia", "ZWE": "zimbabwe", "MOZ": "mozambique", "MDG": "madagascar", "TUN": "tunisia", "MAR": "morocco", "DZA": "algeria",
    "BFA": "burkina faso", "MWI": "malawi", "RWA": "rwanda", "NER": "niger", "AGO": "angola", "SDN": "sudan", "SSD": "south sudan", "LBY": "libya", "GAB": "gabon",
    "COG": "congo", "COD": "congo drc", "BEN": "benin", "TGO": "togo", "NAM": "namibia", "BWA": "botswana", "TCD": "chad", "BDI": "burundi", "SLE": "sierra leone",
    "LBR": "liberia", "GIN": "guinea", "MUS": "mauritius", "GMB": "gambia", "SOM": "somalia", "CIV": "cote d'ivoire",
    # Asia
    "CHN": "china", "JPN": "japan", "KOR": "south korea", "IND": "india", "THA": "thailand", "VNM": "vietnam", "BGD": "bangladesh", "IDN": "indonesia", "KHM": "cambodia",
    "TWN": "taiwan", "MNG": "mongolia", "PAK": "pakistan", "IRN": "iran", "ISR": "israel", "MYS": "malaysia", "PHL": "philippines", "NPL": "nepal", "KAZ": "kazakhstan",
    "LAO": "laos", "MMR": "myanmar", "HKG": "china", "SGP": "singapore", "LKA": "sri lanka", "SAU": "saudi arabia", "TUR": "turkey", "IRQ": "iraq", "LBN": "lebanon",
    "JOR": "jordan", "KWT": "kuwait", "OMN": "oman", "QAT": "qatar", "ARE": "united arab emirates", "BTN": "bhutan", "UZB": "uzbekistan", "AFG": "afghanistan", "YEM": "yemen", "SYR": "syria",
    # Europe
    "DEU": "germany", "FRA": "france", "ITA": "italy", "ESP": "spain", "GBR": "united kingdom", "NLD": "netherlands", "POL": "poland", "BEL": "belgium", "SWE": "sweden",
    "DNK": "denmark", "HUN": "hungary", "CZE": "czech republic", "AUT": "austria", "CHE": "switzerland", "PRT": "portugal", "GRC": "greece", "FIN": "finland", "NOR": "norway",
    "IRL": "ireland", "SRB": "serbia", "RUS": "russia", "ROU": "romania", "BGR": "bulgaria", "HRV": "croatia", "SVN": "slovenia", "SVK": "slovakia", "UKR": "ukraine",
    "BLR": "belarus", "LTU": "lithuania", "LVA": "latvia", "EST": "estonia", "LUX": "luxembourg", "ISL": "iceland", "CYP": "cyprus", "MLT": "malta", "BIH": "bosnia and herzegovina",
    "ALB": "albania", "MDA": "moldova", "MKD": "north macedonia", "MNE": "montenegro",
    # Americas
    "USA": "usa", "CAN": "canada", "MEX": "mexico", "GTM": "guatemala", "CUB": "cuba", "CRI": "costa rica", "PAN": "panama", "HND": "honduras", "SLV": "el salvador",
    "NIC": "nicaragua", "HTI": "haiti", "DOM": "dominican republic", "PRI": "puerto rico", "JAM": "jamaica", "TTO": "trinidad and tobago", "BLZ": "belize",
    "BRA": "brazil", "ARG": "argentina", "PER": "peru", "CHL": "chile", "COL": "colombia", "VEN": "venezuela", "BOL": "bolivia", "ECU": "ecuador", "URY": "uruguay",
    "PRY": "paraguay", "GUF": "french guiana", "SUR": "suriname", "GUY": "guyana",
    # Oceania
    "AUS": "australia", "NZL": "new zealand", "PNG": "papua new guinea", "FJI": "fiji", "NCL": "new caledonia", "WSM": "samoa", "TON": "tonga", "VUT": "vanuatu",
}
