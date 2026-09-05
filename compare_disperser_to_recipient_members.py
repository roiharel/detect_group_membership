from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BINS = Path("outputs/disperser_finescale_integration/disperser_2min_contact_rows.parquet")
EVENTS = Path("outputs/canonical_group_merge_scale_log_scatter/canonical_disperser_events.csv")
MEMBERSHIP = Path(r"C:\Users\rharel\Documents\New project\outputs\canonical_robust_hourly_membership_local_2h_support\canonical_hourly_membership.parquet")
GPS = Path(r"C:\Users\rharel\Documents\New project\outputs\network_v1_cleaned_gps_v1.parquet")
OUT = Path("outputs/disperser_vs_recipient_members")
RADII = (2.0, 5.0)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare disperser proximity integration with recipient-group residents in the same bins.")
    p.add_argument("--bins-file", type=Path, default=BINS); p.add_argument("--events-file", type=Path, default=EVENTS)
    p.add_argument("--membership-file", type=Path, default=MEMBERSHIP); p.add_argument("--gps-file", type=Path, default=GPS)
    p.add_argument("--output-dir", type=Path, default=OUT); p.add_argument("--bootstrap-replicates", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20260808); return p.parse_args()


def distances(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    r=6_371_000.0; lo=np.radians(lon); la=np.radians(lat)
    dlo=lo[:,None]-lo[None,:]; dla=la[:,None]-la[None,:]
    a=np.sin(dla/2)**2+np.cos(la[:,None])*np.cos(la[None,:])*np.sin(dlo/2)**2
    return r*2*np.arctan2(np.sqrt(a),np.sqrt(np.maximum(0,1-a)))


def bootstrap(cells: pd.DataFrame, n: int, seed: int) -> pd.DataFrame:
    ids=cells.event_id.drop_duplicates().to_numpy(); rng=np.random.default_rng(seed); out=[]
    for rep in range(n):
        chosen=rng.choice(ids,len(ids),replace=True)
        s=pd.concat([cells[cells.event_id.eq(x)] for x in chosen],ignore_index=True)
        g=s.groupby(["radius_m","progress_bin"],observed=True)[["focal_contact_probability","resident_contact_probability"]].mean().reset_index()
        g["replicate"]=rep; g["gap_focal_minus_resident"]=g.focal_contact_probability-g.resident_contact_probability; out.append(g)
    return pd.concat(out,ignore_index=True)


def plot(summary: pd.DataFrame, ci: pd.DataFrame, animal: pd.DataFrame, output: Path) -> None:
    colors={"focal":"#7b3294","resident":"#777777"}; fig,axes=plt.subplots(1,2,figsize=(15,6))
    for radius,style in [(2.0,"--"),(5.0,"-")]:
        f=summary[summary.radius_m.eq(radius)].sort_values('progress_bin'); c=ci[ci.radius_m.eq(radius)].set_index('progress_bin').reindex(f.progress_bin)
        x=(f.progress_bin+.5)*10
        for who,label in [("focal","Disperser"),("resident","Recipient members")]:
            y=f[f"{who}_contact_probability"]; axes[0].plot(x,y,style,color=colors[who],marker='o',label=f"{label}, {int(radius)} m")
            axes[0].fill_between(x,c[f"{who}_low"],c[f"{who}_high"],color=colors[who],alpha=.08)
    axes[0].set_xlabel('Progress through dispersal event (%)');axes[0].set_ylabel('Contact probability');axes[0].set_title('Change through event time');axes[0].legend(fontsize=8,frameon=False)
    a=animal[animal.radius_m.eq(5)].sort_values('gap_focal_minus_resident'); labels=[f"{r.animal_id}: {r.origin_group} → {r.recipient_group}" for r in a.itertuples()]
    axes[1].barh(labels,a.gap_focal_minus_resident,color=np.where(a.gap_focal_minus_resident>=0,'#7b3294','#777777'))
    axes[1].axvline(0,color='#444444',lw=1);axes[1].set_xlabel('Disperser minus resident contact probability');axes[1].set_title('Relative integration at 5 m')
    fig.suptitle('Disperser integration versus established recipient-group members\nSame two-minute bins; event-bootstrap 95% CIs')
    fig.tight_layout();fig.savefig(output,dpi=220,bbox_inches='tight');plt.close(fig)


def main() -> None:
    args=parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    focal=pd.read_parquet(args.bins_file); focal=focal[focal.radius_m.isin(RADII)].copy();focal.bin_2min=pd.to_datetime(focal.bin_2min);focal.hour=pd.to_datetime(focal.hour)
    base=focal[["event_id","animal_id","origin_group","recipient_group","hour","bin_2min"]].drop_duplicates()
    event_hours=base[["event_id","animal_id","recipient_group","hour"]].drop_duplicates()
    membership=pd.read_parquet(args.membership_file,columns=["animal_id","window_start","dynamic_social_unit"]);membership.window_start=pd.to_datetime(membership.window_start)
    member_hours=membership.merge(event_hours,left_on='window_start',right_on='hour',suffixes=('','_focal'),how='inner')
    member_hours=member_hours[member_hours.dynamic_social_unit.eq(member_hours.recipient_group)&~member_hours.animal_id.eq(member_hours.animal_id_focal)]
    member_hours=member_hours[["event_id","hour","animal_id"]].drop_duplicates()
    animals=sorted(member_hours.animal_id.astype(str).unique()); start=base.hour.min().tz_localize('UTC');end=(base.hour.max()+pd.Timedelta(hours=1)).tz_localize('UTC')
    gps=pd.read_parquet(args.gps_file,columns=['animal_id','timestamp','location.long','location.lat'],filters=[('animal_id','in',animals),('timestamp','>=',start),('timestamp','<',end)])
    gps.timestamp=pd.to_datetime(gps.timestamp,utc=True).dt.tz_localize(None);gps['hour']=gps.timestamp.dt.floor('h');gps['bin_2min']=gps.timestamp.dt.floor('2min')
    gps=gps.rename(columns={'location.long':'longitude','location.lat':'latitude'}).groupby(['hour','bin_2min','animal_id'],observed=True,as_index=False).agg(longitude=('longitude','median'),latitude=('latitude','median'))
    residents=gps.merge(member_hours,on=['hour','animal_id'],how='inner').merge(base[['event_id','bin_2min']].drop_duplicates(),on=['event_id','bin_2min'],how='inner')
    resident_rows=[]
    for (event_id,bin_time),g in residents.groupby(['event_id','bin_2min'],observed=True):
        n=len(g); d=distances(g.longitude.to_numpy(),g.latitude.to_numpy()) if n>1 else np.empty((n,n))
        for radius in RADII:
            if n<2: rate=np.nan
            else:
                adjacency=(d<=radius)&~np.eye(n,dtype=bool);rate=float(adjacency.any(axis=1).mean())
            resident_rows.append({'event_id':event_id,'bin_2min':bin_time,'radius_m':radius,'resident_n':n,'resident_contact_fraction':rate})
    resident=pd.DataFrame(resident_rows)
    paired=focal.merge(resident,on=['event_id','bin_2min','radius_m'],how='left')
    events=pd.read_csv(args.events_file,parse_dates=['start_time','end_time'])[['event_id','start_time','end_time']];paired=paired.merge(events,on='event_id',how='left')
    dur=(paired.end_time-paired.start_time).dt.total_seconds().clip(lower=3600);paired['progress_bin']=np.floor(((paired.bin_2min-paired.start_time).dt.total_seconds()/dur).clip(0,.999999)*10).astype(int)
    cells=paired.groupby(['event_id','animal_id','origin_group','recipient_group','radius_m','progress_bin'],observed=True).agg(
        focal_contact_probability=('recipient_contact','mean'),resident_contact_probability=('resident_contact_fraction','mean'),bins=('bin_2min','size')).reset_index()
    valid=cells.dropna(subset=['resident_contact_probability']);summary=valid.groupby(['radius_m','progress_bin'],observed=True).agg(
        events=('event_id','nunique'),focal_contact_probability=('focal_contact_probability','mean'),resident_contact_probability=('resident_contact_probability','mean')).reset_index()
    draws=bootstrap(valid,args.bootstrap_replicates,args.seed);q=draws.groupby(['radius_m','progress_bin']).quantile([.025,.975]).unstack()
    ci=pd.DataFrame({'focal_low':q[('focal_contact_probability',.025)],'focal_high':q[('focal_contact_probability',.975)],'resident_low':q[('resident_contact_probability',.025)],'resident_high':q[('resident_contact_probability',.975)]}).reset_index()
    animal=valid.groupby(['animal_id','origin_group','recipient_group','radius_m'],observed=True).agg(events=('event_id','nunique'),focal_contact_probability=('focal_contact_probability','mean'),resident_contact_probability=('resident_contact_probability','mean')).reset_index();animal['gap_focal_minus_resident']=animal.focal_contact_probability-animal.resident_contact_probability
    paired.to_parquet(args.output_dir/'disperser_recipient_member_bin_comparison.parquet',index=False);cells.to_csv(args.output_dir/'disperser_recipient_member_event_progress_cells.csv',index=False)
    summary.to_csv(args.output_dir/'disperser_recipient_member_progress_summary.csv',index=False);animal.to_csv(args.output_dir/'disperser_recipient_member_animal_summary.csv',index=False)
    plot(summary,ci,animal,args.output_dir/'disperser_vs_recipient_members.png')
    meta={'events_with_resident_comparison':int(valid.event_id.nunique()),'animals':int(valid.animal_id.nunique()),'radii_m':RADII,'resident_benchmark':'mean proportion of simultaneously observed recipient residents within radius of another recipient resident','bootstrap_cluster_unit':'dispersal event'}
    (args.output_dir/'disperser_vs_recipient_members_metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8');print(json.dumps(meta,indent=2));print(animal[animal.radius_m.eq(5)].to_string(index=False))


if __name__=='__main__':main()
