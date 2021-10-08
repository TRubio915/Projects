using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ItemToPick : MonoBehaviour
{
  public Item thisItem;
  public InventoryScript playerInventory;

    public void Awake()
    {
        
    }


    private void OnTriggerEnter2D(Collider2D other)
    {
        if (other.gameObject.CompareTag("Player"))
          {
            AddNewItem();
            Destroy(gameObject);
        }
    }

    private void AddNewItem()
    {
        if (!playerInventory.itemList.Contains(thisItem))
        {

            playerInventory.itemList.Add(thisItem);
            //thisItem.itemAmount += 1;
            //InventoryManager.CreateNewItem(thisItem);
        }
        else
        {
            //thisItem.itemAmount += 1;
        }
        InventoryManager.RefreshItem();
    }
}
