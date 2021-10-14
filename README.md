# 内容

本コードはスローンデジタルスカイサーベイ(SDSS)のLuminous Red Galaxy (LRG)カタログからハローカタログを再構築するコードです。

SDSS/LRGは広い赤方偏移にわたって均質な銀河サンプルが得られるため宇宙論的な情報を引き出すのに広く用いられています。
一方、理論側ではN体シミュレーションを通してダークマターが集積して重力的な平衡状態に達したダークマターハロー(以下ハローとよぶ)の研究が進み、
ハローモデルとよばれる現象論的なモデルが確立しており、ハローとLRGをどのように対応づけるかが課題となります。

一般に銀河はハロー中心付近に位置する大質量のセントラル銀河とハロー内をランダムに運動するサテライト銀河に大別することができます。
LRGの多くはセントラル銀河であると考えられますが、サテライト銀河も一部含まれており、それを取り除くことができればハローとの対応関係が非常に理解しやすくなります。

今回、Counts-in-Cylinder法を使って近接したLRGをグループ化し、グループ内で最も明るい銀河のみを選択することでハローに対応するカタログを作成しました。
ハローカタログを解析したところ、サテライト銀河によって引き起こされる非線形な赤方偏移変形(Fingers-of-God効果)の影響が見られず、
サテライト銀河の多くが取り除かれていることを確認しました。

# References

- Where are the Luminous Red Galaxies (LRGs)? Using correlation measurements and lensing to relate LRGs to dark matter halos  
Chiaki Hikage, Rachel Mandelbaum, Masahiro Takada, David N. Spergel  
Mon. Not. Roy. Astron. Soc., Vol.435, Issue 3 (2013) 2345-2370

- Impacts of satellite galaxies on the redshift-space distortions  
Chiaki Hikage, Kazuhiro Yamamoto  
J. Cosmol. Astropart. Phys. Vol 2013, Issue 08, id.19
