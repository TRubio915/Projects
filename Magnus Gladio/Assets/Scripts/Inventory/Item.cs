using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[CreateAssetMenu(fileName ="New Item",menuName ="Inventory/New Item" )]
public class Item : ScriptableObject
{
  public string itemName;
  public Sprite itemImage;
  //public int itemAmount;
  public string itemType;
  [TextArea]
  public string itemInfo;

  public bool consumable;
   
}
